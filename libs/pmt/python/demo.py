import time
import pmt


@pmt.measure("nvidia")
def my_kernel1():
    time.sleep(1)


def my_kernel2():
    platform = "nvidia"
    pm = pmt.create(platform, 0)
    start = pm.read()
    time.sleep(1)
    end = pm.read()
    return {
        "platform": platform,
        "joules": format(pmt.joules(start, end), ".3f"),
        "seconds": format(pmt.seconds(start, end), ".3f"),
        "watt": format(pmt.watts(start, end), ".3f"),
    }


if __name__ == "__main__":
    print(my_kernel1())
    print(my_kernel2())
