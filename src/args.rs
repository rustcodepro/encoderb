use clap::{Parser, Subcommand};
#[derive(Debug, Parser)]
#[command(
    name = "encoderb",
    version = "1.0",
    about = "sparse kmer approach for classification of variant.
       ************************************************
       Gaurav Sablok,
       Email: codeprog@icloud.com
      ************************************************"
)]
pub struct CommandParse {
    /// subcommands for the specific actions
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Classify the variant sequences
    Classifier {
        /// inputfile for the sequences containing the variants
        filepathinput: String,
        /// input file for the sequences on which you want to see association
        predictinput: String,
        /// sparsekmer
        kmer: String,
    },
}
