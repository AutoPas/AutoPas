from .pmt import create, PMT
import re


def is_tty_device(device_id):
    return isinstance(device_id, str) and bool(re.search(r"^/dev/tty.*", device_id))


def joules(start, end):
    return pmt.PMT.joules(start, end)


def seconds(start, end):
    return pmt.PMT.seconds(start, end)


def watts(start, end):
    return pmt.PMT.watts(start, end)


def measure(platform, device_id=None):
    device_id = str(device_id) if device_id else ""

    def decorator(func):
        def wrapper(*args, **kwargs):
            pm = pmt.create(platform, device_id)

            start = pm.read()
            func_results = func(*args, **kwargs)
            end = pm.read()

            results = []
            if func_results is not None:
                results.append(func_results)
            pmt_results = {
                "platform": platform,
                "joules": format(joules(start, end), ".3f"),
                "seconds": format(seconds(start, end), ".3f"),
                "watt": format(watts(start, end), ".3f"),
            }
            results.append(pmt_results)

            return results

        return wrapper

    return decorator


def dump(platform, **kwargs):
    filename = None
    filename_arg = "filename"
    if filename_arg in kwargs:
        filename = kwargs[filename_arg]

    device_id = 0
    device_id_arg = "device_id"
    if device_id_arg in kwargs:
        device_id = kwargs[device_id_arg]

    def decorator(func):
        def wrapper(*args, **kwargs):
            pm = pm.create(platform, device_id)
            pm.startDump(filename)
            result = func(*args, **kwargs)
            pm.stopDump()
            return result

        return wrapper

    return decorator
