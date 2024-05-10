import sys
import time
import traceback
import typing

FUNC_NAMES = {'count_tads_change_intensity': ['Searching changes in the intensity of TADs...',
                                              'The number of TADs with a changed intensity: '],
              'main_split_merge_detection': ['Searching splits and merges in TADs...',
                                             'The number of splits|merges: ']
              }


def wrapper_print(func: typing.Callable) -> typing.Callable:
    def wrapper(*args, **kwargs) -> typing.NoReturn:
        sys.stderr.write(f'{FUNC_NAMES[func.__name__][0]}\r')
        sys.stderr.flush()
        try:
            result = func(*args, **kwargs)
            sys.stderr.write(f'Completed successfully!                       \r')
            sys.stderr.flush()
            time.sleep(2)
            sys.stdout.write(f'{FUNC_NAMES[func.__name__][1]}{result}\n')
            sys.stdout.flush()
            time.sleep(2)
        except Exception as e:
            sys.stderr.write(traceback.format_exc())

    return wrapper
