OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg reg_measure[5];
rz(-2.3440517111951347) q[0];
sx q[0];
rz(5.687761401231555) q[0];
sx q[0];
rz(6.595618377006447) q[0];
rz(-1.9662714480176942) q[1];
sx q[1];
rz(3.207299140163332) q[1];
sx q[1];
rz(10.010063414417731) q[1];
cx q[1],q[0];
rz(-0.5428378684580615) q[0];
sx q[0];
rz(4.49518340855523) q[0];
sx q[0];
rz(9.081695748704515) q[0];
rz(-pi/2) q[1];
sx q[1];
rz(3.6167479051307616) q[1];
sx q[1];
rz(5*pi/2) q[1];
cx q[1],q[0];
rz(2.0057018099940827) q[0];
sx q[0];
rz(5.0820281579160325) q[0];
sx q[0];
rz(6.449475399405039) q[0];
rz(-pi/2) q[1];
sx q[1];
rz(3.493962505840196) q[1];
sx q[1];
rz(2*pi) q[1];
cx q[1],q[0];
rz(-3.101690044733331) q[0];
sx q[0];
rz(4.851261323006206) q[0];
sx q[0];
rz(12.229293130560993) q[0];
rz(-1.7576205456558416) q[0];
sx q[0];
rz(4.7402242234419) q[0];
sx q[0];
rz(11.122058396540483) q[0];
rz(-2.586670204557496) q[1];
sx q[1];
rz(4.909935008941431) q[1];
sx q[1];
rz(11.781875383888464) q[1];
rz(-0.6231158050386796) q[1];
sx q[1];
rz(4.7055961812252205) q[1];
sx q[1];
rz(10.116579956181607) q[1];
rz(0.060881542882380124) q[2];
sx q[2];
rz(3.464082235537097) q[2];
sx q[2];
rz(11.47601759128422) q[2];
rz(-2.045766209041801) q[3];
sx q[3];
rz(3.67664711061566) q[3];
sx q[3];
rz(10.197225527059944) q[3];
rz(0.8951429396595607) q[4];
sx q[4];
rz(4.363765387114951) q[4];
sx q[4];
rz(12.367853176697016) q[4];
cx q[4],q[3];
rz(-0.9476683033986415) q[3];
sx q[3];
rz(4.378455964777682) q[3];
sx q[3];
rz(9.193451480977885) q[3];
rz(-pi/2) q[4];
sx q[4];
rz(3.905370718085974) q[4];
sx q[4];
rz(5*pi/2) q[4];
cx q[4],q[3];
rz(2.0057018099940827) q[3];
sx q[3];
rz(5.0820281579160325) q[3];
sx q[3];
rz(6.449475399405039) q[3];
rz(-pi/2) q[4];
sx q[4];
rz(3.822629388635926) q[4];
sx q[4];
rz(2*pi) q[4];
cx q[4],q[3];
rz(1.6200415432162565) q[3];
sx q[3];
rz(4.01468042002978) q[3];
sx q[3];
rz(8.547850736787739) q[3];
rz(1.4325292567005539) q[3];
sx q[3];
rz(4.839001083632585) q[3];
sx q[3];
rz(8.80406079078967) q[3];
cx q[3],q[0];
rz(-0.6309311124018606) q[0];
sx q[0];
rz(4.465515225509809) q[0];
sx q[0];
rz(9.101968582708198) q[0];
rz(-pi/2) q[3];
sx q[3];
rz(3.7031066075990324) q[3];
sx q[3];
rz(5*pi/2) q[3];
cx q[3],q[0];
rz(2.0057018099940827) q[0];
sx q[0];
rz(5.0820281579160325) q[0];
sx q[0];
rz(6.449475399405039) q[0];
rz(pi/2) q[3];
sx q[3];
rz(3.4831126567276414) q[3];
sx q[3];
rz(3*pi) q[3];
cx q[3],q[0];
rz(-0.8456301596876452) q[0];
sx q[0];
rz(5.2525294581245445) q[0];
sx q[0];
rz(12.538608830202115) q[0];
rz(-0.5734203513494807) q[0];
sx q[0];
rz(4.350126605771845) q[0];
sx q[0];
rz(8.088260961986744) q[0];
rz(-1.9486197017440992) q[3];
sx q[3];
rz(5.09007048727638) q[3];
sx q[3];
rz(8.235912297287108) q[3];
rz(0.8129041982048246) q[3];
sx q[3];
rz(5.336632061533088) q[3];
sx q[3];
rz(10.95938096893338) q[3];
cx q[3],q[1];
rz(-0.8686592675606679) q[1];
sx q[1];
rz(4.397103421951131) q[1];
sx q[1];
rz(9.16823958031011) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(3.99243627097477) q[3];
sx q[3];
rz(5*pi/2) q[3];
cx q[3],q[1];
rz(2.0057018099940827) q[1];
sx q[1];
rz(5.0820281579160325) q[1];
sx q[1];
rz(6.449475399405039) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(3.686173806529009) q[3];
sx q[3];
rz(2*pi) q[3];
cx q[3],q[1];
rz(-2.7597738652199904) q[1];
sx q[1];
rz(3.832318466604609) q[1];
sx q[1];
rz(12.292153012951625) q[1];
rz(2.0942159885524108) q[1];
sx q[1];
rz(3.5501264510778023) q[1];
sx q[1];
rz(9.086697807303345) q[1];
rz(1.4289461974163604) q[3];
sx q[3];
rz(3.4131358804222174) q[3];
sx q[3];
rz(10.39384698325474) q[3];
rz(2.830405286394438) q[3];
sx q[3];
rz(3.834358640803329) q[3];
sx q[3];
rz(6.996316740829839) q[3];
rz(2.826603734267353) q[4];
sx q[4];
rz(5.191224551750153) q[4];
sx q[4];
rz(11.951560288160087) q[4];
rz(-1.2406344015938275) q[4];
sx q[4];
rz(4.785127321360154) q[4];
sx q[4];
rz(7.445201549403353) q[4];
cx q[2],q[4];
rz(-pi/2) q[2];
sx q[2];
rz(4.27844969832182) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.6774405316944794) q[4];
sx q[4];
rz(4.450746982587107) q[4];
sx q[4];
rz(9.113669262053694) q[4];
cx q[2],q[4];
rz(-pi/2) q[2];
sx q[2];
rz(3.35960013483742) q[2];
sx q[2];
rz(2*pi) q[2];
rz(2.0057018099940827) q[4];
sx q[4];
rz(5.0820281579160325) q[4];
sx q[4];
rz(6.449475399405039) q[4];
cx q[2],q[4];
rz(-0.9807160046969625) q[2];
sx q[2];
rz(5.735544348610964) q[2];
sx q[2];
rz(12.475189791232093) q[2];
rz(1.7539125803435294) q[2];
sx q[2];
rz(3.801375750732422) q[2];
sx q[2];
rz(11.396230826731939) q[2];
rz(-0.30632906923761327) q[4];
sx q[4];
rz(5.402465243429531) q[4];
sx q[4];
rz(8.201988520201228) q[4];
rz(2.9070400607744196) q[4];
sx q[4];
rz(4.039722550881752) q[4];
sx q[4];
rz(7.325950323180919) q[4];
cx q[4],q[2];
rz(-0.7348173479016769) q[2];
sx q[2];
rz(4.433428732866076) q[2];
sx q[2];
rz(9.128994592623613) q[2];
rz(-pi/2) q[4];
sx q[4];
rz(3.706243919428896) q[4];
sx q[4];
rz(5*pi/2) q[4];
cx q[4],q[2];
rz(2.0057018099940827) q[2];
sx q[2];
rz(5.0820281579160325) q[2];
sx q[2];
rz(6.449475399405039) q[2];
rz(pi/2) q[4];
sx q[4];
rz(3.3140518317948615) q[4];
sx q[4];
rz(3*pi) q[4];
cx q[4],q[2];
rz(2.4157337622281396) q[2];
sx q[2];
rz(4.709944866007869) q[2];
sx q[2];
rz(7.308261122850187) q[2];
rz(-1.211258381794331) q[2];
sx q[2];
rz(3.785790190020495) q[2];
sx q[2];
rz(10.792775782875607) q[2];
cx q[2],q[3];
rz(-pi/2) q[2];
sx q[2];
rz(3.726999092799904) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.8517380182323908) q[3];
sx q[3];
rz(4.401371290267834) q[3];
sx q[3];
rz(9.163026720938815) q[3];
cx q[2],q[3];
rz(-pi/2) q[2];
sx q[2];
rz(3.438073368997294) q[2];
sx q[2];
rz(2*pi) q[2];
rz(2.0057018099940827) q[3];
sx q[3];
rz(5.0820281579160325) q[3];
sx q[3];
rz(6.449475399405039) q[3];
cx q[2],q[3];
rz(1.5399442817063065) q[2];
sx q[2];
rz(3.9317918077224747) q[2];
sx q[2];
rz(12.401462121246555) q[2];
rz(2.8556791655661193) q[2];
sx q[2];
rz(3.6365591886710504) q[2];
sx q[2];
rz(8.690564499738702) q[2];
rz(-3.0921544845294475) q[3];
sx q[3];
rz(4.821660684541987) q[3];
sx q[3];
rz(6.73775547769181) q[3];
rz(2.435004846284909) q[3];
sx q[3];
rz(4.673529690148774) q[3];
sx q[3];
rz(7.432863290666097) q[3];
cx q[3],q[2];
rz(-0.6188273631964236) q[2];
sx q[2];
rz(4.469462461300634) q[2];
sx q[2];
rz(9.099033874983396) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(3.730395701507332) q[3];
sx q[3];
rz(5*pi/2) q[3];
cx q[3],q[2];
rz(2.0057018099940827) q[2];
sx q[2];
rz(5.0820281579160325) q[2];
sx q[2];
rz(6.449475399405039) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(3.1659315474358554) q[3];
sx q[3];
rz(2*pi) q[3];
cx q[3],q[2];
rz(0.4505277683183384) q[2];
sx q[2];
rz(4.897676241225634) q[2];
sx q[2];
rz(6.390245484055367) q[2];
rz(2.959454790727518) q[3];
sx q[3];
rz(5.7345886876338845) q[3];
sx q[3];
rz(8.1662398061698) q[3];
rz(-0.6292721071010052) q[4];
sx q[4];
rz(5.01178525270341) q[4];
sx q[4];
rz(8.31469429762215) q[4];
rz(-0.8497506569029527) q[4];
sx q[4];
rz(5.818698757255057) q[4];
sx q[4];
rz(6.846344944702436) q[4];
cx q[4],q[0];
rz(-0.6996786992713329) q[0];
sx q[0];
rz(4.443914633520314) q[0];
sx q[0];
rz(9.119495150337665) q[0];
rz(-pi/2) q[4];
sx q[4];
rz(3.8795679705012276) q[4];
sx q[4];
rz(5*pi/2) q[4];
cx q[4],q[0];
rz(2.0057018099940827) q[0];
sx q[0];
rz(5.0820281579160325) q[0];
sx q[0];
rz(6.449475399405039) q[0];
rz(pi/2) q[4];
sx q[4];
rz(3.7624701216759586) q[4];
sx q[4];
rz(3*pi) q[4];
cx q[4],q[0];
rz(-0.2975479032472452) q[0];
sx q[0];
rz(5.61079442213053) q[0];
sx q[0];
rz(11.545735099509788) q[0];
rz(-0.7633079216740515) q[4];
sx q[4];
rz(4.190783729566778) q[4];
sx q[4];
rz(7.975022171107442) q[4];
rz(-1.2224992395485137) q[4];
sx q[4];
rz(4.948825301816922) q[4];
sx q[4];
rz(7.0395812666265005) q[4];
cx q[1],q[4];
rz(-pi/2) q[1];
sx q[1];
rz(4.275066143086883) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.6510694231525385) q[4];
sx q[4];
rz(4.459042155022315) q[4];
sx q[4];
rz(9.106953216747332) q[4];
cx q[1],q[4];
rz(-pi/2) q[1];
sx q[1];
rz(3.357263790115109) q[1];
sx q[1];
rz(2*pi) q[1];
rz(2.0057018099940827) q[4];
sx q[4];
rz(5.0820281579160325) q[4];
sx q[4];
rz(6.449475399405039) q[4];
cx q[1],q[4];
rz(1.1538321707063588) q[1];
sx q[1];
rz(4.605835623762798) q[1];
sx q[1];
rz(6.368411677865275) q[1];
rz(-0.5092454218621709) q[4];
sx q[4];
rz(5.769647154455248) q[4];
sx q[4];
rz(11.380805274888836) q[4];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> reg_measure[0];
measure q[1] -> reg_measure[1];
measure q[2] -> reg_measure[2];
measure q[3] -> reg_measure[3];
measure q[4] -> reg_measure[4];
