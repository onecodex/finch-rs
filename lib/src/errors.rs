use std::result::Result as StdResult;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum FinchError {
    #[error("failed to load/read/write file: {0:?}")]
    Io(#[from] std::io::Error),
    #[error("capnproto error: {0:?}")]
    Capnproto(#[from] capnp::Error),
    #[error("failed to parse the fasta/fastq file: {0}")]
    Needletail(#[from] needletail::errors::ParseError),
    #[error("failed to parse as integer")]
    IntError(#[from] core::num::ParseIntError),
    #[error("failed to parse as float")]
    FloatError(#[from] core::num::ParseFloatError),
    #[error("enum value not found in schema")]
    SchemaError(#[from] capnp::NotInSchema),
    #[error("json error: {0:?}")]
    Json(#[from] serde_json::Error),
    #[error("Finch error: {0}")]
    Message(String),
}

pub type FinchResult<T> = StdResult<T, FinchError>;

// TODO: remove the macro magic if possible when done with moving off failure

#[doc(hidden)]
#[macro_export]
macro_rules! bail {
    ($e:expr) => {
        return Err($crate::errors::FinchError::Message($e.to_owned()));
    };
    ($fmt:expr, $($arg:tt)*) => {
        return Err($crate::errors::FinchError::Message(format!($fmt, $($arg)*)))
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! format_err {
    ($($arg:tt)*) => { $crate::errors::FinchError::Message(format!($($arg)*)) }
}
