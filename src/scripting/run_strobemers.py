import os

def run_strobemer(
    protocol,
    order,
    strobe_length,
    w_min,
    w_max,
    experiment_name,
    experiment_version
):
    command = (
        "cargo run --bin strobemer_comparison -- "
        f"references.fasta "
        f"query.fasta "
        f"-p {protocol} "
        f"-o {order} "
        f"-l {strobe_length} "
        f"--w-min {w_min} "
        f"--w-max {w_max} "
        f"-e {experiment_name} "
        f"""-v {experiment_version} """
    )
    print(f"Executing: {command}")
    os.system(command)

def strobemer_experiment1():
    run_strobemer(
        protocol="rand",
        order=2,
        strobe_length=5,
        w_min=1, w_max=50,
        experiment_name="experiment1",
        experiment_version=""
    )


def strobemer_experiment1_1():
    run_strobemer(
        protocol="rand",
        order=2,
        strobe_length=10,
        w_min=1, w_max=15,
        experiment_name="experiment1",
        experiment_version="_1"
    )


def strobemer_experiment1_2():
    run_strobemer(
        protocol="rand",
        order=3,
        strobe_length=10,
        w_min=1, w_max=50,
        experiment_name="experiment1",
        experiment_version="_2"
    )

def strobemer_experiment1_3():
    run_strobemer(
        protocol="min",
        order=2,
        strobe_length=5,
        w_min=0, w_max=40,
        experiment_name="experiment1",
        experiment_version="_3"
    )




def strobemer_experiment2_1():
    run_strobemer(
        protocol="rand",
        order=2,
        strobe_length=15,
        w_min=0, w_max=30,
        experiment_name="experiment2",
        experiment_version="_1"
    )


def strobemer_experiment2_2():
    run_strobemer(
        protocol="rand",
        order=3,
        strobe_length=10,
        w_min=0, w_max=25,
        experiment_name="experiment2",
        experiment_version="_2"
    )

def strobemer_experiment2_3():
    run_strobemer(
        protocol="rand",
        order=2,
        strobe_length=15,
        w_min=0, w_max=30,
        experiment_name="experiment2",
        experiment_version="_3"
    )

def strobemer_experiment2_4():
    run_strobemer(
        protocol="min",
        order=2,
        strobe_length=15,
        w_min=0, w_max=30,
        experiment_name="experiment2",
        experiment_version="_4"
    )


def strobemer_experiment3_1():
    run_strobemer(
        protocol="rand",
        order=2,
        strobe_length=15,
        w_min=0, w_max=30,
        experiment_name="experiment3",
        experiment_version="_1"
    )

strobemer_experiment1_1()