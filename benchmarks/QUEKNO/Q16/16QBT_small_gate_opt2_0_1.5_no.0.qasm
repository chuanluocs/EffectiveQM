OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
h q[7];
cx q[11],q[7];
cx q[3],q[11];
h q[2];
h q[3];
cx q[13],q[9];
cx q[13],q[2];
h q[7];
h q[2];
cx q[11],q[7];
h q[13];
h q[7];
h q[7];
