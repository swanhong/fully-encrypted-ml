use clap::Parser;
use std::fmt;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Argument {
    /// target algorithm
    #[arg(short, long, default_value = "protocol")]
    pub target_algorithm: String,

    /// bit_len for the group
    /// (default = 100)
    /// (options: 10, 100, 3072)
    #[arg(short, long, default_value_t = 100)]
    pub bit_len: u64,

    /// input dimension for the protocol
    /// (default = 5)
    #[arg(long, default_value_t = 5)]
    pub dim0: usize,

    /// intermediate dimension for the protocol
    /// (default = 4)
    #[arg(long, default_value_t = 4)]
    pub dim1: usize,

    /// output dimension for the protocol
    /// (default = 3)
    #[arg(long, default_value_t = 3)]
    pub dim2: usize,

    /// message & function bound
    /// (default = 10)
    #[arg(long, default_value_t = 10)]
    pub bound: usize,
    
    /// number of trys to check runtime and correctness of the protocol
    #[arg(short, long, default_value_t = 1)]
    pub n_try: usize,

    /// scaling facor for ml algorithm
    /// (default = 1073741824)
    #[arg(short, long, default_value_t = 1073741824)]
    pub scale: i64,

    /// decomposition number
    /// (default = 4)
    #[arg(long, default_value_t = 4)]
    pub n_decomp: usize,
    
    /// (optional) data_ml if you want to run ml example (options: breast, iris)
    /// (default = breast)
    #[arg(short, long, default_value = "breast")]
    pub data_ml: Option<String>,
}

// impl Argument {
//     pub fn new() -> Argument {
//         Argument {
//             target_algorithm: "protocol".to_string(),
//             dim0: 5,
//             dim1: 4,
//             dim2: 3,
//             n_try: 1,
//             data_ml: Some("breast".to_string()),
//             bound: 10,
//             bit_len: 100,
//             scale: 1073741824,
//             n_decomp: 4,
//         }
//     }
// }

impl fmt::Display for Argument {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "=== pQuant Parameters ===\n\
             target_algorithm   : {}\n\
             bit_len            : {}\n\
             n_decomp           : {}\n",
            self.target_algorithm, 
            self.bit_len,
            self.n_decomp,
        )?;
        if self.target_algorithm == "protocol" {
            write!(
                f,
                "dim0               : {}\n\
                 dim1               : {}\n\
                 dim2               : {}\n\
                 bound              : {}\n\
                 n_try              : {}\n",
                self.dim0,
                self.dim1,
                self.dim2,
                self.bound,
                self.n_try,
            )
        } else if self.target_algorithm == "ml" {
            write!(
                f,
                "scale              : {}\n\
                 data_ml            : {}\n",
                self.scale,
                self.data_ml.as_ref().map(|s| s.clone()).unwrap_or_default(),
            )
        } else {
            write!(
                f,
                "=== pQuant Parameters ===\n\
                 target_algorithm   : {}\n\
                 bit_len            : {}\n\
                 dim0               : {}\n\
                 dim1               : {}\n\
                 dim2               : {}\n\
                 bound              : {}\n\
                 n_try              : {}\n\
                 data_ml            : {}\n\
                 scale              : {}\n\
                 n_decomp           : {}",
                self.target_algorithm,
                self.bit_len,
                self.dim0,
                self.dim1,
                self.dim2,
                self.bound,
                self.n_try,
                self.data_ml.as_ref().map(|s| s.clone()).unwrap_or_default(),
                self.scale,
                self.n_decomp,
            )
        }
    }
}
