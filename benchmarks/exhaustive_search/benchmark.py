from subprocess import Popen, PIPE
import time

for device in ["cpu", "gpu"]:
    for nodes in range(18,24):
        input_file = f"mvp{nodes}.qubo"
        output_file = f"/tmp/result{nodes}.csv"
        run = f'../../build/src/one-solver-exhaustive --input examples/{input_file} --output {output_file} --device={device}'.split(" ")

        start = time.time()
        p = Popen(run, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate()
        rc = p.returncode
        elapsed_time = (time.time() - start)

        if rc==0:
            device_name = output.split(b'\n')[-2].decode("ascii")
            print(f'nodes={nodes}, device_name="{device_name}", elapsed_time={elapsed_time}')