def _print_progress_bar(
    iteration: int,
    total: int,
    *,
    prefix: str = "",
    suffix: str = "",
    decimals: int = 1,
    length: int = 100,
    fill: str = "*",
    print_end: str = "\r",
):
    r"""Call in a loop to create terminal progress bar.

    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(
        100 if (total == 0) else 100 * (iteration / float(total))
    )
    filled_length = int(length if (total == 0) else length * iteration // total)
    bar = fill * filled_length + "-" * (length - filled_length)
    print(f"\r{prefix} |{bar}| {percent}% {suffix}", end=print_end)
    # Print New Line on Complete
    if iteration == total:
        print()


class ProgressBar:
    def __init__(self, message: str, max_n: int = 100) -> None:
        self.n = max_n
        self.i = 0
        self.message = message
        _print_progress_bar(self.i, self.n, prefix=self.message, length=50)

    def next(self) -> None:
        self.i = self.i + 1
        _print_progress_bar(self.i, self.n, prefix=self.message, length=50)

    def finish(self) -> None:
        _print_progress_bar(self.n, self.n, prefix=self.message, length=50)
