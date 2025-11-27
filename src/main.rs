use std::fs::FileTimes;

use self::encoderb::Sequences;
use self::encoderb::VariantSeq;
use crate::args::CommandParse;
use crate::args::Commands;
use clap::Parser;
use smartcore::metrics::*;
use smartcore::svm::svc::SVC;
use smartcore::svm::svc::SVCParameters;
mod args;
mod encoderb;

/*
 Gaurav Sablok
 codeprog@icloud.com
*/

fn main() {
    let argparse = CommandParse::parse();
    match &argparse.command {
        Commands::Classifier {
            filepathinput,
            predictinput,
            kmer,
        } => {
            let sequencepath = Sequences {
                path: filepathinput.to_string(),
            };
            let classlabels = sequencepath.sequences().unwrap();
            let prediction = Sequences {
                path: predictinput.to_string(),
            };
            let predictionlabel = prediction.sequences().unwrap();
            let train_dataset = VariantSeq::generatesparse(filepathinput, classlabels.1, kmer);
            let test_dataset = VariantSeq::new(classlabels.0, classlabels.1, kmer);
            let train_data = train_dataset.dataset();
            let test_data = test_dataset.dataset();
            let smartcore_run = SVC::fit(&train_data, &test_data, SVCParameters::with_c(self, 1));
            let model = smartcore_run.fit(&train_data).expect("Failed to train SVM");
            let train_predictions = model.predict(&train_data);
            let train_accuracy = accuracy(&test_data, &train_predictions).round();
            let test_predictions = model.predict(&test_data);
            let test_accuracy = accuracy(&test_dataset, &test_predictions).round();
        }
    }
}
