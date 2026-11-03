//! Child-process jobs of `polytool atlas`: every item of a seeding or growth
//! run is computed by `polytool atlas job` in a child process with a hard
//! time limit (the parent kills it on expiry), or, for batches that share a
//! parent library, in a watched thread of one child with a per-item time
//! limit.  The parent integrates the results in a deterministic order.

use polytool::atlas::grow::{forward_derivations, ForwardTask};
use polytool::atlas::seed::{
    basis_map_theorem_derivation, certify_definition_equality, lift_context, lift_derivations,
    root_derivations, AnalyzedDerivation, LiftContext,
};
use polytool::atlas::{AtlasCertificate, AtlasSequenceSpec, AtlasStore};
use polytool::recurrence::RecurrenceJson;
use serde::{Deserialize, Serialize};
use std::io::{Read, Write};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{mpsc, Arc, Mutex};
use std::time::{Duration, Instant};

/// A request to `polytool atlas job` (read from stdin as JSON).
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "job", rename_all = "snake_case")]
pub(crate) enum JobRequest {
    /// Root derivations of one source (and its basis-map theorem).
    Roots { spec: AtlasSequenceSpec },
    /// Lift derivations of several sources, sharing the libraries built
    /// from the atlas in `directory`.
    Lifts {
        directory: PathBuf,
        specs: Vec<AtlasSequenceSpec>,
        item_timeout_seconds: u64,
    },
    /// An all-index equality certificate.
    Equality {
        existing: RecurrenceJson,
        new: RecurrenceJson,
    },
    /// Forward tasks.
    Forward {
        tasks: Vec<ForwardTask>,
        item_timeout_seconds: u64,
    },
}

/// The result of one item of a job.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct JobItem {
    /// `ok`, `timeout`, or `error: ...`.
    pub(crate) status: String,
    /// Wall-clock seconds.
    pub(crate) seconds: f64,
    /// The derivations found.
    pub(crate) derivations: Vec<AnalyzedDerivation>,
    /// For an equality job: the certificate, when certified.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) equality: Option<AtlasCertificate>,
}

impl JobItem {
    pub(crate) fn failed(status: impl Into<String>, seconds: f64) -> Self {
        Self {
            status: status.into(),
            seconds,
            derivations: Vec::new(),
            equality: None,
        }
    }
}

/// The response of a job (one item per source or task).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct JobResponse {
    pub(crate) items: Vec<JobItem>,
}

/// Run `f` in a thread with a deadline; `None` on expiry (the thread is
/// left to finish in the background of this child process, which exits
/// when its job is done).
fn with_deadline<T: Send + 'static>(
    timeout: Duration,
    f: impl FnOnce() -> T + Send + 'static,
) -> Option<T> {
    let (sender, receiver) = mpsc::channel();
    let spawned = std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(move || {
            let _ = sender.send(f());
        });
    if spawned.is_err() {
        return None;
    }
    receiver.recv_timeout(timeout).ok()
}

fn item_from(result: Result<Vec<AnalyzedDerivation>, String>, start: Instant) -> JobItem {
    let seconds = start.elapsed().as_secs_f64();
    match result {
        Ok(derivations) => JobItem {
            status: "ok".into(),
            seconds,
            derivations,
            equality: None,
        },
        Err(error) => JobItem::failed(format!("error: {error}"), seconds),
    }
}

/// Execute a job in this process (the `polytool atlas job` command).
pub(crate) fn execute_job(request: JobRequest) -> JobResponse {
    match request {
        JobRequest::Roots { spec } => {
            let start = Instant::now();
            let result = root_derivations(&spec)
                .map_err(|e| e.to_string())
                .map(|mut out| {
                    if !out.is_empty() {
                        if let Ok(Some(map)) = basis_map_theorem_derivation(&spec) {
                            out.push(map);
                        }
                    }
                    out
                });
            JobResponse {
                items: vec![item_from(result, start)],
            }
        }
        JobRequest::Equality { existing, new } => {
            let start = Instant::now();
            let item = match certify_definition_equality(&existing, &new) {
                Ok(certificate) => JobItem {
                    status: "ok".into(),
                    seconds: start.elapsed().as_secs_f64(),
                    derivations: Vec::new(),
                    equality: Some(certificate),
                },
                Err(error) => {
                    JobItem::failed(format!("error: {error}"), start.elapsed().as_secs_f64())
                }
            };
            JobResponse { items: vec![item] }
        }
        JobRequest::Lifts {
            directory,
            specs,
            item_timeout_seconds,
        } => {
            let context: Result<Arc<LiftContext>, String> = AtlasStore::open(&directory)
                .and_then(|store| lift_context(&store))
                .map(Arc::new)
                .map_err(|e| e.to_string());
            let items = specs
                .into_iter()
                .map(|spec| {
                    let start = Instant::now();
                    let context = match &context {
                        Ok(context) => context.clone(),
                        Err(error) => return JobItem::failed(format!("error: {error}"), 0.0),
                    };
                    match with_deadline(Duration::from_secs(item_timeout_seconds), move || {
                        lift_derivations(&spec, &context).map_err(|e| e.to_string())
                    }) {
                        Some(result) => item_from(result, start),
                        None => JobItem::failed("timeout", start.elapsed().as_secs_f64()),
                    }
                })
                .collect();
            JobResponse { items }
        }
        JobRequest::Forward {
            tasks,
            item_timeout_seconds,
        } => {
            let items = tasks
                .into_iter()
                .map(|task| {
                    let start = Instant::now();
                    match with_deadline(Duration::from_secs(item_timeout_seconds), move || {
                        forward_derivations(&task).map_err(|e| e.to_string())
                    }) {
                        Some(result) => item_from(result, start),
                        None => JobItem::failed("timeout", start.elapsed().as_secs_f64()),
                    }
                })
                .collect();
            JobResponse { items }
        }
    }
}

/// `polytool atlas job`: read a request from stdin, write the response to
/// stdout.
pub(crate) fn cmd_atlas_job() {
    exit_when_orphaned();
    let mut input = String::new();
    if std::io::stdin().read_to_string(&mut input).is_err() {
        eprintln!("atlas job: cannot read the request");
        std::process::exit(2);
    }
    let request: JobRequest = match serde_json::from_str(&input) {
        Ok(request) => request,
        Err(error) => {
            eprintln!("atlas job: invalid request: {error}");
            std::process::exit(2);
        }
    };
    let response = execute_job(request);
    let text = serde_json::to_string(&response).unwrap_or_else(|_| "{\"items\":[]}".into());
    let mut stdout = std::io::stdout();
    let _ = stdout.write_all(text.as_bytes());
    let _ = stdout.flush();
    std::process::exit(0);
}

/// Exit this child when its parent process is gone (a parent killed by an
/// outer time limit must not leave its jobs running).
fn exit_when_orphaned() {
    #[cfg(unix)]
    {
        let parent = std::os::unix::process::parent_id();
        std::thread::spawn(move || loop {
            std::thread::sleep(Duration::from_millis(500));
            if std::os::unix::process::parent_id() != parent {
                std::process::exit(3);
            }
        });
    }
}

/// Run one job in a child process, killed after `timeout`.
pub(crate) fn run_child(request: &JobRequest, timeout: Duration) -> Result<JobResponse, String> {
    let executable = std::env::current_exe().map_err(|e| e.to_string())?;
    let mut child = Command::new(executable)
        .args(["atlas", "job"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| e.to_string())?;
    let text = serde_json::to_string(request).map_err(|e| e.to_string())?;
    let mut stdin = child.stdin.take().ok_or("no stdin")?;
    let writer = std::thread::spawn(move || {
        let _ = stdin.write_all(text.as_bytes());
    });
    let mut stdout = child.stdout.take().ok_or("no stdout")?;
    let reader = std::thread::spawn(move || {
        let mut output = Vec::new();
        let _ = stdout.read_to_end(&mut output);
        output
    });
    let deadline = Instant::now() + timeout;
    let status = loop {
        match child.try_wait() {
            Ok(Some(status)) => break Some(status),
            Ok(None) if Instant::now() >= deadline => {
                let _ = child.kill();
                let _ = child.wait();
                break None;
            }
            Ok(None) => std::thread::sleep(Duration::from_millis(20)),
            Err(error) => return Err(error.to_string()),
        }
    };
    let _ = writer.join();
    let output = reader.join().unwrap_or_default();
    match status {
        None => Err("timeout".into()),
        Some(status) if !status.success() => Err(format!("child failed: {status}")),
        Some(_) => serde_json::from_slice(&output).map_err(|e| format!("bad response: {e}")),
    }
}

/// Run the jobs with `parallel` children at a time; results in input order.
pub(crate) fn run_pool(
    requests: Vec<JobRequest>,
    parallel: usize,
    timeout_of: impl Fn(&JobRequest) -> Duration + Sync,
    label: &str,
) -> Vec<Result<JobResponse, String>> {
    let total = requests.len();
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    let results: Mutex<Vec<Option<Result<JobResponse, String>>>> =
        Mutex::new((0..total).map(|_| None).collect());
    let started = Instant::now();
    std::thread::scope(|scope| {
        for _ in 0..parallel.max(1) {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= total {
                    break;
                }
                let result = run_child(&requests[i], timeout_of(&requests[i]));
                results.lock().expect("results")[i] = Some(result);
                let finished = done.fetch_add(1, Ordering::SeqCst) + 1;
                if finished % 25 == 0 || finished == total {
                    eprintln!(
                        "  {label}: {finished}/{total} jobs, {:.0} s",
                        started.elapsed().as_secs_f64()
                    );
                }
            });
        }
    });
    results
        .into_inner()
        .expect("results")
        .into_iter()
        .map(|r| r.unwrap_or_else(|| Err("not run".into())))
        .collect()
}
