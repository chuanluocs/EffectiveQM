OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
cx q[11],q[10];
h q[0];
h q[0];
h q[5];
cx q[0],q[15];
h q[0];
h q[5];
cx q[9],q[5];
cx q[9],q[5];
h q[0];