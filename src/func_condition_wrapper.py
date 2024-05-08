import sys
import traceback
import logging

FUNC_NAMES = {'count_tads_change_intensity': ['Searching changes in the intensity of TADs...',
                                              'The number of TADs with a changed intensity: '],
              'main_split_merge_detection': ['Searching splits and merges in TADs...',
                                             'The number of splits|merges: ']
              }


def wrapper_print(func):
    def wrapper(*args, **kwargs):
        sys.stdout.write(f'{FUNC_NAMES[func.__name__][0]}\n'); sys.stdout.flush()
        try:
            result = func(*args, **kwargs)
            sys.stdout.write(f'Completed successfully! {FUNC_NAMES[func.__name__][1]}{result}\n'); sys.stdout.flush()
        except Exception as e:
            sys.stderr.write(traceback.format_exc())

    return wrapper
