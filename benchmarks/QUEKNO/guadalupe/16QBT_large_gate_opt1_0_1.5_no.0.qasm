OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
cx q[2],q[14];
cx q[11],q[1];
cx q[1],q[7];
h q[6];
h q[1];
cx q[15],q[11];
h q[15];
h q[6];
h q[6];
h q[14];
h q[14];
cx q[15],q[11];
h q[12];
cx q[11],q[1];
h q[13];
cx q[12],q[13];
h q[2];
h q[11];
cx q[14],q[12];
h q[7];
h q[12];
cx q[11],q[6];
h q[13];