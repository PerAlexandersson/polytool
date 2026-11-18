//! Persistent content-addressed reports and bounded JSON views.
//! Only validated report ids are accepted; callers cannot supply read paths.

use polytool::digest::sha256_hex;
use serde_json::{json, Map, Value};
use std::fs::{self, File, FileTimes, OpenOptions};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::time::{Duration, SystemTime};

const DEFAULT_MAX_BYTES: u64 = 512 * 1024 * 1024;
const DEFAULT_MAX_AGE: u64 = 14 * 24 * 60 * 60;
const MAX_REPORT_BYTES: u64 = 64 * 1024 * 1024;
static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Report storage is opt-in through a directory flag or environment variable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReportStoreConfig {
    pub directory: Option<PathBuf>,
    pub max_bytes: u64,
    pub max_age: Duration,
}

impl Default for ReportStoreConfig {
    fn default() -> Self {
        Self {
            directory: None,
            max_bytes: DEFAULT_MAX_BYTES,
            max_age: Duration::from_secs(DEFAULT_MAX_AGE),
        }
    }
}

impl ReportStoreConfig {
    pub fn from_env() -> Result<Self, String> {
        fn number(name: &str, default: u64) -> Result<u64, String> {
            match std::env::var(name) {
                Ok(v) => v
                    .parse()
                    .map_err(|_| format!("{name} must be a nonnegative integer")),
                Err(std::env::VarError::NotPresent) => Ok(default),
                Err(_) => Err(format!("{name} must be Unicode")),
            }
        }
        Ok(Self {
            directory: std::env::var_os("POLYTOOL_MCP_REPORT_DIR")
                .filter(|p| !p.is_empty())
                .map(PathBuf::from),
            max_bytes: number("POLYTOOL_MCP_REPORT_MAX_BYTES", DEFAULT_MAX_BYTES)?,
            max_age: Duration::from_secs(number(
                "POLYTOOL_MCP_REPORT_MAX_AGE_SECONDS",
                DEFAULT_MAX_AGE,
            )?),
        })
    }
}

#[derive(Debug)]
struct StoreInner {
    config: ReportStoreConfig,
    writes: AtomicU64,
}

#[derive(Debug, Clone)]
pub(crate) struct ReportStore(Arc<StoreInner>);

fn error(message: &str) -> io::Error {
    io::Error::other(message)
}

fn valid_id(id: &str) -> bool {
    id.len() == 34
        && id.starts_with("r-")
        && id[2..]
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

fn report_id(bytes: &[u8]) -> String {
    format!("r-{}", &sha256_hex(bytes)[..32])
}

// Refuse symlink directories as well as symlink report files. Never follow a
// caller-controlled link when creating, reading, touching or evicting reports.
fn regular_directory(path: &Path) -> io::Result<()> {
    for ancestor in path.ancestors() {
        match fs::symlink_metadata(ancestor) {
            Ok(m) if m.file_type().is_symlink() || !m.is_dir() => {
                return Err(error(
                    "report directory must contain only regular directories",
                ))
            }
            Ok(_) => {}
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => return Err(e),
        }
    }
    Ok(())
}

impl ReportStore {
    pub(crate) fn new(config: ReportStoreConfig) -> Self {
        let store = Self(Arc::new(StoreInner {
            config,
            writes: AtomicU64::new(0),
        }));
        // Failure is reported with the next computed summary, never as a
        // computation failure or an unexpected startup failure.
        let _ = store.evict();
        store
    }

    fn directory(&self) -> io::Result<&Path> {
        let directory = self
            .0
            .config
            .directory
            .as_deref()
            .ok_or_else(|| error("report store is disabled; configure a report directory"))?;
        if directory.to_str().is_none() {
            return Err(error(
                "report directory must be representable as a JSON path",
            ));
        }
        regular_directory(directory)?;
        Ok(directory)
    }

    pub(crate) fn store(&self, bytes: &[u8]) -> Option<Value> {
        self.0.config.directory.as_ref()?;
        Some(match self.write(bytes) {
            Ok(value) => value,
            Err(_) => {
                json!({"stored":false,"reason":"report store write failed; check directory permissions and storage limits"})
            }
        })
    }

    fn write(&self, bytes: &[u8]) -> io::Result<Value> {
        let directory = self.directory()?;
        fs::create_dir_all(directory)?;
        if bytes.len() as u64 > MAX_REPORT_BYTES || bytes.len() as u64 > self.0.config.max_bytes {
            return Err(error("report exceeds store capacity"));
        }
        let id = report_id(bytes);
        let path = directory.join(format!("{id}.json"));
        match fs::symlink_metadata(&path) {
            Ok(metadata) => {
                if !metadata.is_file() || self.read(&id)? != bytes {
                    return Err(error("existing report is corrupt or not a regular file"));
                }
            }
            Err(e) if e.kind() == io::ErrorKind::NotFound => {
                let temporary = directory.join(format!(
                    "{id}.json.tmp-{}-{}",
                    std::process::id(),
                    TEMP_COUNTER.fetch_add(1, Ordering::Relaxed)
                ));
                let result = (|| {
                    let mut file = OpenOptions::new()
                        .write(true)
                        .create_new(true)
                        .open(&temporary)?;
                    file.write_all(bytes)?;
                    file.sync_all()?;
                    // Another writer may have finished the same content. Its
                    // regular file is left untouched; otherwise rename atomically.
                    match fs::symlink_metadata(&path) {
                        Ok(m) if m.is_file() => {
                            fs::remove_file(&temporary)?;
                        }
                        Ok(_) => return Err(error("report destination is not a regular file")),
                        Err(e) if e.kind() == io::ErrorKind::NotFound => {
                            fs::rename(&temporary, &path)?
                        }
                        Err(e) => return Err(e),
                    }
                    Ok(())
                })();
                if result.is_err() {
                    let _ = fs::remove_file(&temporary);
                }
                result?;
            }
            Err(e) => return Err(e),
        }
        if self.0.writes.fetch_add(1, Ordering::Relaxed) % 64 == 63 {
            self.evict()?;
        }
        Ok(
            json!({"id":id,"uri":format!("polytool://reports/{id}"),"path":path,"bytes":bytes.len()}),
        )
    }

    pub(crate) fn read(&self, id: &str) -> io::Result<Vec<u8>> {
        if !valid_id(id) {
            return Err(error(
                "invalid report id; expected r- followed by 32 lowercase hexadecimal digits",
            ));
        }
        let directory = self.directory()?;
        let path = directory.join(format!("{id}.json"));
        let missing = || error("report evicted or unknown; re-run (results are deterministic)");
        let metadata = fs::symlink_metadata(&path).map_err(|e| {
            if e.kind() == io::ErrorKind::NotFound {
                missing()
            } else {
                e
            }
        })?;
        if !metadata.is_file() || metadata.len() > MAX_REPORT_BYTES {
            return Err(error("report is not a regular file or exceeds read limit"));
        }
        let file = File::open(path)?;
        let mut bytes = vec![];
        (&file).take(MAX_REPORT_BYTES + 1).read_to_end(&mut bytes)?;
        if bytes.len() as u64 > MAX_REPORT_BYTES || report_id(&bytes) != id {
            return Err(error(
                "report hash mismatch; corrupted report must be recomputed",
            ));
        }
        file.set_times(FileTimes::new().set_modified(SystemTime::now()))?;
        Ok(bytes)
    }

    fn evict(&self) -> io::Result<()> {
        if self.0.config.directory.is_none() {
            return Ok(());
        }
        let directory = self.directory()?;
        if !directory.exists() {
            return Ok(());
        }
        let now = SystemTime::now();
        let mut files = vec![];
        for entry in fs::read_dir(directory)? {
            let entry = entry?;
            let name = entry.file_name();
            let Some(id) = name.to_str().and_then(|s| s.strip_suffix(".json")) else {
                continue;
            };
            if !valid_id(id) || !entry.file_type()?.is_file() {
                continue;
            }
            let m = entry.metadata()?;
            files.push((m.modified()?, entry.path(), m.len()));
        }
        files.sort();
        let mut total: u64 = files.iter().map(|f| f.2).sum();
        for (mtime, path, bytes) in files {
            if total > self.0.config.max_bytes
                || now.duration_since(mtime).unwrap_or_default() > self.0.config.max_age
            {
                match fs::remove_file(path) {
                    Ok(()) => total = total.saturating_sub(bytes),
                    Err(e) if e.kind() == io::ErrorKind::NotFound => {
                        total = total.saturating_sub(bytes)
                    }
                    Err(e) => return Err(e),
                }
            }
        }
        Ok(())
    }
}

pub(crate) fn validate_selection(
    select: Option<&[String]>,
    max_bytes: Option<usize>,
) -> Result<usize, String> {
    let limit = max_bytes.unwrap_or(8192);
    if !(128..=65536).contains(&limit) {
        return Err("max_bytes must be between 128 and 65536".into());
    }
    if let Some(pointers) = select {
        if pointers.len() > 64 {
            return Err("select accepts at most 64 pointers".into());
        }
        for p in pointers {
            pointer_segments(p)?;
        }
    }
    Ok(limit)
}

fn pointer_segments(pointer: &str) -> Result<Vec<String>, String> {
    if pointer.is_empty() {
        return Ok(vec![]);
    }
    if !pointer.starts_with('/') || pointer.len() > 1024 {
        return Err("select requires RFC 6901 pointers of at most 1024 bytes".into());
    }
    let mut wildcards = 0;
    let mut segments = vec![];
    for segment in pointer[1..].split('/') {
        let mut decoded = String::new();
        let mut chars = segment.chars();
        while let Some(c) = chars.next() {
            if c == '~' {
                decoded.push(match chars.next() {
                    Some('0') => '~',
                    Some('1') => '/',
                    _ => return Err("invalid RFC 6901 escape".into()),
                });
            } else {
                decoded.push(c);
            }
        }
        if decoded == "*" {
            wildcards += 1;
        }
        segments.push(decoded);
    }
    if wildcards > 1 {
        return Err("select permits at most one wildcard segment per pointer".into());
    }
    Ok(segments)
}

fn escaped(segment: &str) -> String {
    segment.replace('~', "~0").replace('/', "~1")
}

pub(crate) fn selected(value: &Value, pointers: &[String], limit: usize) -> (Value, bool) {
    fn visit(
        value: &Value,
        segments: &[String],
        path: String,
        output: &mut Map<String, Value>,
        limit: usize,
        truncated: &mut bool,
    ) {
        if *truncated {
            return;
        }
        let Some((first, rest)) = segments.split_first() else {
            output.insert(path.clone(), value.clone());
            if serde_json::to_vec(&output).unwrap().len() > limit {
                output.remove(&path);
                *truncated = true;
            }
            return;
        };
        if first == "*" {
            match value {
                Value::Array(items) => {
                    for (i, v) in items.iter().enumerate() {
                        visit(v, rest, format!("{path}/{i}"), output, limit, truncated);
                        if *truncated {
                            break;
                        }
                    }
                }
                Value::Object(items) => {
                    for (k, v) in items {
                        visit(
                            v,
                            rest,
                            format!("{path}/{}", escaped(k)),
                            output,
                            limit,
                            truncated,
                        );
                        if *truncated {
                            break;
                        }
                    }
                }
                _ => {}
            }
        } else if let Some(child) = value.pointer(&format!("/{}", escaped(first))) {
            visit(
                child,
                rest,
                format!("{path}/{}", escaped(first)),
                output,
                limit,
                truncated,
            );
        }
    }
    let mut output = Map::new();
    let mut truncated = false;
    for p in pointers {
        // Validated by the request entry point before running computations.
        let segments = pointer_segments(p).expect("validated pointer");
        visit(
            value,
            &segments,
            String::new(),
            &mut output,
            limit,
            &mut truncated,
        );
    }
    (Value::Object(output), truncated)
}

pub(crate) fn outline(value: &Value) -> (Value, bool) {
    fn view(v: &Value, depth: usize) -> Value {
        match v {
            Value::Array(a) => json!({"length":a.len()}),
            Value::Object(o) if depth == 0 => json!({"keys":o.keys().collect::<Vec<_>>()}),
            Value::Object(o) => Value::Object(
                o.iter()
                    .map(|(k, v)| (k.clone(), view(v, depth - 1)))
                    .collect(),
            ),
            _ => v.clone(),
        }
    }
    let mut output = Map::new();
    let mut truncated = false;
    if let Some(object) = value.as_object() {
        for (key, v) in object {
            output.insert(key.clone(), view(v, 1));
            if serde_json::to_vec(&output).unwrap().len() > 1800 {
                output.remove(key);
                truncated = true;
            }
        }
    }
    (Value::Object(output), truncated)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub(crate) struct TempDirectory(pub PathBuf);
    impl TempDirectory {
        pub(crate) fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "polytool-report-test-{}-{}",
                std::process::id(),
                TEMP_COUNTER.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }
        pub(crate) fn config(&self) -> ReportStoreConfig {
            ReportStoreConfig {
                directory: Some(self.0.clone()),
                ..Default::default()
            }
        }
    }
    impl Drop for TempDirectory {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.0).unwrap();
        }
    }

    #[test]
    fn stable_handles_survive_restart_and_concurrent_identical_writes() {
        let temp = TempDirectory::new();
        let store = ReportStore::new(temp.config());
        let bytes = b"{\"schema\":\"test\",\"verified\":true}";
        let first = store.store(bytes).unwrap();
        let second = ReportStore::new(temp.config());
        assert_eq!(second.store(bytes).unwrap(), first);
        let threads: Vec<_> = (0..8)
            .map(|_| {
                let store = ReportStore::new(temp.config());
                std::thread::spawn(move || store.store(bytes).unwrap())
            })
            .collect();
        for t in threads {
            assert_eq!(t.join().unwrap(), first);
        }
        assert_eq!(fs::read_dir(&temp.0).unwrap().count(), 1);
        assert_eq!(second.read(first["id"].as_str().unwrap()).unwrap(), bytes);
        assert!(fs::symlink_metadata(first["path"].as_str().unwrap())
            .unwrap()
            .is_file());
    }

    #[test]
    fn rejects_traversal_missing_and_corrupted_reports_and_touches_reads() {
        let temp = TempDirectory::new();
        let store = ReportStore::new(temp.config());
        for id in [
            "../escape",
            "r-0000000000000000000000000000000/",
            "r-AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        ] {
            assert!(store
                .read(id)
                .unwrap_err()
                .to_string()
                .contains("invalid report id"));
        }
        assert!(store
            .read("r-00000000000000000000000000000000")
            .unwrap_err()
            .to_string()
            .contains("evicted or unknown"));
        let h = store.store(b"{}").unwrap();
        let path = h["path"].as_str().unwrap();
        let id = h["id"].as_str().unwrap();
        let old = SystemTime::now() - Duration::from_secs(1000);
        File::open(path)
            .unwrap()
            .set_times(FileTimes::new().set_modified(old))
            .unwrap();
        store.read(id).unwrap();
        assert!(fs::metadata(path).unwrap().modified().unwrap() > old);
        fs::write(path, b"bad").unwrap();
        assert!(store
            .read(id)
            .unwrap_err()
            .to_string()
            .contains("hash mismatch"));
        assert_eq!(store.store(b"{}").unwrap()["stored"], false);
    }

    #[test]
    fn eviction_uses_age_and_lru_size_at_startup_and_every_64_writes() {
        let temp = TempDirectory::new();
        let store = ReportStore::new(temp.config());
        let older = store.store(b"{\"n\":1}").unwrap();
        let newer = store.store(b"{\"n\":2}").unwrap();
        let old_path = older["path"].as_str().unwrap();
        File::open(old_path)
            .unwrap()
            .set_times(FileTimes::new().set_modified(SystemTime::now() - Duration::from_secs(100)))
            .unwrap();
        let config = ReportStoreConfig {
            max_bytes: 7,
            ..temp.config()
        };
        let small = ReportStore::new(config);
        assert!(!Path::new(old_path).exists());
        assert!(small.read(newer["id"].as_str().unwrap()).is_ok());
        let aged = ReportStore::new(ReportStoreConfig {
            max_age: Duration::ZERO,
            ..temp.config()
        });
        aged.evict().unwrap();
        assert_eq!(fs::read_dir(&temp.0).unwrap().count(), 0);
        let limited = ReportStore::new(ReportStoreConfig {
            max_bytes: 14,
            ..temp.config()
        });
        for i in 0..64 {
            limited.store(format!("{i:07}").as_bytes()).unwrap();
        }
        assert!(fs::read_dir(&temp.0).unwrap().count() <= 2);
    }

    #[test]
    fn disabled_and_failed_writes_leave_computation_available() {
        let disabled = ReportStore::new(ReportStoreConfig::default());
        assert!(disabled.store(b"{}").is_none());
        assert!(disabled.read("r-00000000000000000000000000000000").is_err());
        let temp = TempDirectory::new();
        let path = temp.0.join("file");
        fs::write(&path, b"unchanged").unwrap();
        let failed = ReportStore::new(ReportStoreConfig {
            directory: Some(path.clone()),
            ..Default::default()
        });
        assert_eq!(failed.store(b"{}").unwrap()["stored"], false);
        assert_eq!(fs::read(path).unwrap(), b"unchanged");
    }

    #[test]
    fn pointers_expand_once_decode_escapes_and_bound_large_values() {
        let v = json!({"a/b":{"~":[{"status":"yes"},{"status":"no"}]},"large":"x".repeat(9000)});
        let pointers = vec!["/a~1b/~0/*/status".into()];
        validate_selection(Some(&pointers), None).unwrap();
        let (s, truncated) = selected(&v, &pointers, 1000);
        assert!(!truncated);
        assert_eq!(s["/a~1b/~0/0/status"], "yes");
        assert_eq!(s["/a~1b/~0/1/status"], "no");
        let (s, truncated) = selected(&v, &["/large".into()], 200);
        assert!(truncated);
        assert_eq!(s, json!({}));
        for p in ["bad", "/bad~2", "/*/*"] {
            assert!(validate_selection(Some(&[p.into()]), None).is_err());
        }
        let (view, truncated) = outline(&v);
        assert!(truncated);
        assert!(serde_json::to_vec(&view).unwrap().len() < 2048);
    }
}
