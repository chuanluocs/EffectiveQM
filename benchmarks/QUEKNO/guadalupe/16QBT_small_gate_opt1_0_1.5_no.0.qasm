OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
h q[12];
cx q[0],q[6];
h q[2];
h q[15];
h q[11];
h q[8];
cx q[0],q[6];
cx q[8],q[10];
h q[3];
h q[12];
h q[2];
cx q[14],q[15];
h q[12];
h q[10];
cx q[14],q[3];
cx q[12],q[2];
cx q[14],q[3];
h q[12];
cx q[6],q[11];
h q[0];
