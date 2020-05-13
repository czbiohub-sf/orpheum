import s3fs
import pandas as pd


def write_s3(df, filename, fmt="csv", **kwargs):
    fs = s3fs.S3FileSystem(anon=False)
    if fmt == "csv":
        # csv is a text format
        with fs.open(filename, "w") as f:
            return df.to_csv(f, **kwargs)
    elif fmt == "parquet":
        # Parquet is a binary format and needs the "b" flag
        with fs.open(filename, "wb") as f:
            return df.to_parquet(f, **kwargs)


def read_aws_s3_ls(filename, **kwargs):
    return pd.read_table(
        filename,
        delim_whitespace=True,
        header=None,
        names=["date", "time", "bytes", "basename"],
        **kwargs
    )


def savefig(fig, filename, **kwargs):
    """Saves figure to s3 if path starts with s3:// otherwise saves locally"""
    if filename.startswith("s3://"):
        fs = s3fs.S3FileSystem(anon=False)
        with fs.open(filename, "wb") as f:
            return fig.savefig(f, **kwargs)
    else:
        return fig.savefig(filename)
