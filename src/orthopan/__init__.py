"""OrthoPan command line entry point."""

import argparse


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(prog="orthopan", description="Scalable pangenome exploration")
    subparsers = parser.add_subparsers(dest="command", required=True)
    from orthopan import importer, server

    importer.add_parser(subparsers)
    server.add_parser(subparsers)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
