#!/usr/local/bin/python
import os
import re
import gzip
import argparse


def validate_file(f):
    if not os.path.exists(f):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(f))
    return f


def t2g_line(transcript, gene):
    return str(transcript) + "\t" + str(gene) + "\n"


def build_t2g_re():
    return re.compile("([A-Z]+[0-9]+)\.\S+\s([A-Z]+[0-9]+)\.\S+")


def get_t2g(line, t2g_re):
    return t2g_re.match(line)


def get_t2g_line(line, t2g_re):
    t2g_id = get_t2g(line, t2g_re)
    return {'transcript_id': t2g_id, 'gene_id': t2g_id}


def write_t2g_line(t2g, line, t2g_re):
    results = get_t2g_line(line, t2g_re)
    if results['transcript_id']:
        t2g.write(
            t2g_line(
                results['transcript_id'].group(1),
                results['gene_id'].group(2)
            )
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="create transcript to genes file from a gtf file."
    )
    parser.add_argument(
        "-f", "--t2g", dest="t2g", required=True, type=validate_file,
        help="t2g file", metavar="FILE"
    )
    args = parser.parse_args()
    t2g_re = build_t2g_re()

    try:
        with gzip.open(args.t2g, "rb") as gtf:
            with open("fix_t2g.txt", "w") as t2g:
                for line in gtf:
                    write_t2g_line(t2g, str(line), t2g_re)
    except gzip.BadGzipFile:
        with open(args.t2g, "r") as gtf:
            with open("fix_t2g.txt", "w") as t2g:
                for line in gtf:
                    write_t2g_line(t2g, str(line), t2g_re)
