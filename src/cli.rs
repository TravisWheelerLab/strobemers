use clap::Parser;

#[derive(Debug, Parser)]
#[command(about, author, version)]
pub struct CommonArgs {
    #[arg(required = true, value_name = "REF FILE", help = "References FASTA file")]
    pub references_file: String,
    #[arg(required = true, value_name = "QUERY FILE", help = "Query FAST(A/Q) file")]
    pub query_file: String,
    //#[arg(short='m', value_name = "METHOD()", help = "E.g. jaccard_similarity")] // unused
    //pub similarity_method: String,
}