OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
cx q[4],q[13];
h q[13];
h q[10];
cx q[10],q[12];
h q[3];
cx q[4],q[13];
h q[15];
h q[12];
cx q[3],q[12];
cx q[15],q[9];
h q[4];
h q[9];
h q[15];