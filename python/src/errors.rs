//! The exception hierarchy of the Python package.
//!
//! Every library error (invalid input, an unsupported recurrence, a failed
//! replay) is raised as [`PolytoolError`], a subclass of `ValueError`, or one
//! of its subclasses.  Wrong Python types raise the usual `TypeError`.

use pyo3::create_exception;
use pyo3::exceptions::PyValueError;

create_exception!(
    polytool,
    PolytoolError,
    PyValueError,
    "Base class of the errors raised by polytool (a ValueError): invalid input, \
     unsupported recurrences and library errors."
);
create_exception!(
    polytool,
    VerificationError,
    PolytoolError,
    "A certificate or report did not replay: its recorded data differ from the \
     exact recomputation."
);
create_exception!(
    polytool,
    BudgetExhaustedError,
    PolytoolError,
    "A recurrence search reached max_candidates or max_work before it finished. \
     The attributes candidates_tried, candidates_considered and work_units \
     record how far it got."
);
