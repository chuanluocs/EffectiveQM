OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
h q[11];
h q[12];
h q[12];
cx q[2],q[11];
h q[5];
h q[9];
h q[4];
cx q[16],q[2];
cx q[8],q[14];
cx q[9],q[14];
h q[17];
cx q[16],q[11];
h q[12];
cx q[2],q[11];
h q[17];
h q[11];
h q[17];
h q[15];
cx q[0],q[12];
cx q[12],q[15];
cx q[17],q[8];
h q[4];
cx q[17],q[5];
h q[9];
h q[14];
cx q[4],q[18];
h q[14];
h q[14];
