function der_st_dx4 = der_stima_dx4(lambda,mu,V)
%DER_STIMA_DX4
%    DER_ST_DX4 = DER_STIMA_DX4(LAMBDA,MU,X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    02-Jun-2015 16:23:12
x1=V(1,1);
y1=V(1,2);
z1=V(1,3);

x2=V(2,1);
y2=V(2,2);
z2=V(2,3);

x3=V(3,1);
y3=V(3,2);
z3=V(3,3);

x4=V(4,1);
y4=V(4,2);
z4=V(4,3);

t2 = x2.*y3;
t3 = x4.*y2;
t4 = x3.*y4;
t52 = x3.*y2;
t53 = x2.*y4;
t54 = x4.*y3;
t5 = t2+t3+t4-t52-t53-t54;
t6 = x2.*z3;
t7 = x4.*z2;
t8 = x3.*z4;
t44 = x3.*z2;
t45 = x2.*z4;
t46 = x4.*z3;
t9 = t6+t7+t8-t44-t45-t46;
t11 = y2.*z3;
t12 = y3.*z2;
t27 = y2.*z4;
t28 = y4.*z2;
t29 = y3.*z4;
t30 = y4.*z3;
t10 = t11-t12-t27+t28+t29-t30;
t13 = x1.*y2.*z3;
t14 = x2.*y3.*z1;
t15 = x3.*y1.*z2;
t16 = x1.*y4.*z2;
t17 = x2.*y1.*z4;
t18 = x4.*y2.*z1;
t19 = x1.*y3.*z4;
t20 = x3.*y4.*z1;
t21 = x4.*y1.*z3;
t22 = x2.*y4.*z3;
t23 = x3.*y2.*z4;
t24 = x4.*y3.*z2;
t31 = x1.*y3.*z2;
t32 = x2.*y1.*z3;
t33 = x3.*y2.*z1;
t34 = x1.*y2.*z4;
t35 = x2.*y4.*z1;
t36 = x4.*y1.*z2;
t37 = x1.*y4.*z3;
t38 = x3.*y1.*z4;
t39 = x4.*y3.*z1;
t40 = x2.*y3.*z4;
t41 = x3.*y4.*z2;
t42 = x4.*y2.*z3;
t25 = t13+t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t24-t31-t32-t33-t34-t35-t36-t37-t38-t39-t40-t41-t42;
t26 = z2-z3;
t43 = 1.0./t25;
t47 = y1.*z2;
t48 = y3.*z1;
t55 = y2.*z1;
t56 = y1.*z3;
t49 = t11-t12+t47+t48-t55-t56;
t50 = 1.0./t25.^2;
t51 = y2-y3;
t57 = x1.*y3;
t58 = x4.*y1;
t72 = x3.*y1;
t73 = x1.*y4;
t59 = t4-t54+t57+t58-t72-t73;
t60 = x1.*z3;
t61 = x4.*z1;
t68 = x3.*z1;
t69 = x1.*z4;
t62 = t8-t46+t60+t61-t68-t69;
t63 = mu.*2.0;
t64 = lambda+t63;
t65 = z1-z3;
t66 = y4.*z1;
t70 = y1.*z4;
t67 = t29-t30-t48+t56+t66-t70;
t71 = y1-y3;
t74 = x2.*y1;
t82 = x1.*y2;
t75 = t3-t53-t58+t73+t74-t82;
t76 = x2.*z1;
t80 = x1.*z2;
t77 = t7-t45-t61+t69+t76-t80;
t78 = z1-z2;
t79 = t27-t28+t47-t55+t66-t70;
t81 = y1-y2;
t83 = t2-t52-t57+t72-t74+t82;
t84 = t6-t44-t60+t68-t76+t80;
t85 = lambda.*t10.*t26;
t86 = mu.*t10.*t26;
t87 = t85+t86;
t88 = t43.*t87.*(1.0./6.0);
t89 = lambda.*t9.*t10;
t90 = mu.*t9.*t10;
t91 = t89+t90;
t92 = t49.*t50.*t91.*(1.0./6.0);
t93 = t88+t92;
t94 = mu.*t5.*t51.*2.0;
t95 = t5.^2;
t96 = mu.*t95;
t97 = t10.^2;
t98 = t9.^2;
t99 = mu.*t51.*t59;
t100 = mu.*t5.*t71;
t101 = mu.*t5.*t59;
t102 = mu.*t5.*t81;
t103 = mu.*t5.*t75;
t104 = mu.*t51.*t83;
t105 = mu.*t5.*t83;
t106 = lambda.*t10.*t51;
t107 = mu.*t10.*t51;
t108 = t106+t107;
t109 = lambda.*t5.*t10;
t110 = mu.*t5.*t10;
t111 = t109+t110;
t112 = t43.*t108.*(-1.0./6.0)-t49.*t50.*t111.*(1.0./6.0);
t113 = lambda.*t5.*t26;
t114 = mu.*t5.*t26;
t115 = lambda.*t9.*t51;
t116 = mu.*t9.*t51;
t117 = t113+t114+t115+t116;
t118 = t43.*t117.*(1.0./6.0);
t119 = lambda.*t5.*t9;
t120 = mu.*t5.*t9;
t121 = t119+t120;
t122 = t49.*t50.*t121.*(1.0./6.0);
t123 = t118+t122;
t124 = mu.*t9.*t26.*2.0;
t125 = mu.*t98;
t126 = mu.*t97;
t127 = mu.*t26.*t62;
t128 = mu.*t9.*t65;
t129 = mu.*t9.*t62;
t130 = mu.*t10.*t67;
t131 = mu.*t9.*t78;
t132 = mu.*t9.*t77;
t133 = mu.*t26.*t84;
t134 = mu.*t9.*t84;
t135 = mu.*t10.*t49;
t136 = t99+t100+t127+t128;
t137 = t43.*t136.*(1.0./6.0);
t138 = t10.*t64.*t67;
t139 = t101+t129+t138;
t140 = t49.*t50.*t139.*(1.0./6.0);
t141 = t137+t140;
t142 = lambda.*t26.*t67;
t143 = mu.*t10.*t65;
t144 = t142+t143;
t145 = lambda.*t9.*t67;
t146 = mu.*t10.*t62;
t147 = t145+t146;
t148 = t43.*t144.*(-1.0./6.0)-t49.*t50.*t147.*(1.0./6.0);
t149 = lambda.*t51.*t67;
t150 = mu.*t10.*t71;
t151 = t149+t150;
t152 = t43.*t151.*(1.0./6.0);
t153 = lambda.*t5.*t67;
t154 = mu.*t10.*t59;
t155 = t153+t154;
t156 = t49.*t50.*t155.*(1.0./6.0);
t157 = t152+t156;
t158 = lambda.*t10.*t65;
t159 = mu.*t26.*t67;
t160 = t158+t159;
t161 = lambda.*t10.*t62;
t162 = mu.*t9.*t67;
t163 = t161+t162;
t164 = t43.*t160.*(-1.0./6.0)-t49.*t50.*t163.*(1.0./6.0);
t165 = t26.*t62.*t64;
t166 = t9.*t64.*t65;
t167 = t99+t100+t165+t166;
t168 = t43.*t167.*(1.0./6.0);
t169 = t9.*t62.*t64;
t170 = t101+t130+t169;
t171 = t49.*t50.*t170.*(1.0./6.0);
t172 = t168+t171;
t173 = lambda.*t5.*t65;
t174 = mu.*t26.*t59;
t175 = lambda.*t51.*t62;
t176 = mu.*t9.*t71;
t177 = t173+t174+t175+t176;
t178 = lambda.*t5.*t62;
t179 = mu.*t9.*t59;
t180 = t178+t179;
t181 = t43.*t177.*(-1.0./6.0)-t49.*t50.*t180.*(1.0./6.0);
t182 = lambda.*t65.*t67;
t183 = mu.*t65.*t67;
t184 = t182+t183;
t185 = t43.*t184.*(1.0./6.0);
t186 = lambda.*t62.*t67;
t187 = mu.*t62.*t67;
t188 = t186+t187;
t189 = t49.*t50.*t188.*(1.0./6.0);
t190 = t185+t189;
t191 = mu.*t59.*t71.*2.0;
t192 = t59.^2;
t193 = mu.*t192;
t194 = t67.^2;
t195 = t62.^2;
t196 = mu.*t59.*t81;
t197 = mu.*t59.*t75;
t198 = mu.*t71.*t83;
t199 = mu.*t59.*t83;
t200 = lambda.*t10.*t71;
t201 = mu.*t51.*t67;
t202 = t200+t201;
t203 = t43.*t202.*(1.0./6.0);
t204 = lambda.*t10.*t59;
t205 = mu.*t5.*t67;
t206 = t204+t205;
t207 = t49.*t50.*t206.*(1.0./6.0);
t208 = t203+t207;
t209 = lambda.*t26.*t59;
t210 = mu.*t5.*t65;
t211 = lambda.*t9.*t71;
t212 = mu.*t51.*t62;
t213 = t209+t210+t211+t212;
t214 = lambda.*t9.*t59;
t215 = mu.*t5.*t62;
t216 = t214+t215;
t217 = t43.*t213.*(-1.0./6.0)-t49.*t50.*t216.*(1.0./6.0);
t218 = t51.*t59.*t64;
t219 = t5.*t64.*t71;
t220 = t127+t128+t218+t219;
t221 = t43.*t220.*(1.0./6.0);
t222 = t5.*t59.*t64;
t223 = t129+t130+t222;
t224 = t49.*t50.*t223.*(1.0./6.0);
t225 = t221+t224;
t226 = lambda.*t67.*t71;
t227 = mu.*t67.*t71;
t228 = t226+t227;
t229 = lambda.*t59.*t67;
t230 = mu.*t59.*t67;
t231 = t229+t230;
t232 = t43.*t228.*(-1.0./6.0)-t49.*t50.*t231.*(1.0./6.0);
t233 = lambda.*t59.*t65;
t234 = mu.*t59.*t65;
t235 = lambda.*t62.*t71;
t236 = mu.*t62.*t71;
t237 = t233+t234+t235+t236;
t238 = t43.*t237.*(1.0./6.0);
t239 = lambda.*t59.*t62;
t240 = mu.*t59.*t62;
t241 = t239+t240;
t242 = t49.*t50.*t241.*(1.0./6.0);
t243 = t238+t242;
t244 = mu.*t62.*t65.*2.0;
t245 = mu.*t195;
t246 = mu.*t194;
t247 = mu.*t62.*t78;
t248 = mu.*t62.*t77;
t249 = mu.*t65.*t84;
t250 = mu.*t62.*t84;
t251 = mu.*t49.*t67;
t296 = mu.*t51.*t75;
t353 = mu.*t26.*t77;
t252 = t102+t131-t296-t353;
t253 = t49.*t50.*(t103+t132-t10.*t64.*t79).*(1.0./6.0);
t254 = t253-t43.*t252.*(1.0./6.0);
t255 = lambda.*t26.*t79;
t256 = mu.*t10.*t78;
t257 = t255+t256;
t258 = t43.*t257.*(1.0./6.0);
t259 = lambda.*t9.*t79;
t260 = t259-mu.*t10.*t77;
t261 = t49.*t50.*t260.*(1.0./6.0);
t262 = t258+t261;
t263 = lambda.*t51.*t79;
t264 = mu.*t10.*t81;
t265 = t263+t264;
t266 = lambda.*t5.*t79;
t267 = t266-mu.*t10.*t75;
t268 = t43.*t265.*(-1.0./6.0)-t49.*t50.*t267.*(1.0./6.0);
t314 = mu.*t71.*t75;
t371 = mu.*t65.*t77;
t269 = t196+t247-t314-t371;
t270 = t43.*t269.*(1.0./6.0);
t271 = t197+t248-t64.*t67.*t79;
t272 = t270-t49.*t50.*t271.*(1.0./6.0);
t273 = lambda.*t65.*t79;
t274 = mu.*t67.*t78;
t275 = t273+t274;
t276 = lambda.*t62.*t79;
t277 = t276-mu.*t67.*t77;
t278 = t43.*t275.*(-1.0./6.0)-t49.*t50.*t277.*(1.0./6.0);
t279 = lambda.*t71.*t79;
t280 = mu.*t67.*t81;
t281 = t279+t280;
t282 = t43.*t281.*(1.0./6.0);
t283 = lambda.*t59.*t79;
t284 = t283-mu.*t67.*t75;
t285 = t49.*t50.*t284.*(1.0./6.0);
t286 = t282+t285;
t287 = t3-t53-t58+t73+t74-t82;
t288 = t7-t45-t61+t69+t76-t80;
t289 = lambda.*t10.*t78;
t290 = mu.*t26.*t79;
t291 = t289+t290;
t292 = t43.*t291.*(1.0./6.0);
t293 = lambda.*t10.*t77;
t294 = t293-mu.*t9.*t79;
t295 = t292-t49.*t50.*t294.*(1.0./6.0);
t297 = t9.*t64.*t78;
t298 = t9.*t64.*t77;
t355 = mu.*t10.*t79;
t299 = t49.*t50.*(t103+t298-t355).*(1.0./6.0);
t300 = lambda.*t5.*t78;
t301 = mu.*t9.*t81;
t302 = t300+t301-lambda.*t51.*t77-mu.*t26.*t75;
t303 = t43.*t302.*(1.0./6.0);
t304 = lambda.*t5.*t77;
t305 = mu.*t9.*t75;
t306 = t304+t305;
t307 = t303-t49.*t50.*t306.*(1.0./6.0);
t308 = lambda.*t67.*t78;
t309 = mu.*t65.*t79;
t310 = t308+t309;
t311 = lambda.*t67.*t77;
t312 = t49.*t50.*(t311-mu.*t62.*t79).*(1.0./6.0);
t313 = t312-t43.*t310.*(1.0./6.0);
t315 = t62.*t64.*t78;
t316 = t62.*t64.*t77;
t373 = mu.*t67.*t79;
t317 = t197+t316-t373;
t318 = lambda.*t59.*t78;
t319 = mu.*t62.*t81;
t320 = t318+t319-lambda.*t71.*t77-mu.*t65.*t75;
t321 = lambda.*t59.*t77;
t322 = mu.*t62.*t75;
t323 = t49.*t50.*(t321+t322).*(1.0./6.0);
t324 = t323-t43.*t320.*(1.0./6.0);
t325 = lambda.*t78.*t79;
t326 = mu.*t78.*t79;
t327 = t325+t326;
t328 = t43.*t327.*(1.0./6.0);
t329 = lambda.*t77.*t79;
t330 = mu.*t77.*t79;
t331 = t329+t330;
t332 = t328-t49.*t50.*t331.*(1.0./6.0);
t333 = mu.*t75.*t81.*2.0;
t334 = t3-t53-t58+t73+t74-t82;
t335 = t79.^2;
t336 = t7-t45-t61+t69+t76-t80;
t337 = mu.*t81.*t83;
t338 = mu.*t75.*t83;
t339 = lambda.*t10.*t81;
t340 = mu.*t51.*t79;
t341 = t339+t340;
t342 = lambda.*t10.*t75;
t343 = t49.*t50.*(t342-mu.*t5.*t79).*(1.0./6.0);
t344 = t343-t43.*t341.*(1.0./6.0);
t345 = mu.*t5.*t78;
t346 = lambda.*t9.*t81;
t347 = t345+t346-lambda.*t26.*t75-mu.*t51.*t77;
t348 = t43.*t347.*(1.0./6.0);
t349 = lambda.*t9.*t75;
t350 = mu.*t5.*t77;
t351 = t349+t350;
t352 = t348-t49.*t50.*t351.*(1.0./6.0);
t354 = t5.*t64.*t81;
t356 = t5.*t64.*t75;
t357 = lambda.*t67.*t81;
t358 = mu.*t71.*t79;
t359 = t357+t358;
t360 = t43.*t359.*(1.0./6.0);
t361 = lambda.*t67.*t75;
t362 = t361-mu.*t59.*t79;
t363 = t360-t49.*t50.*t362.*(1.0./6.0);
t364 = lambda.*t65.*t75;
t365 = mu.*t71.*t77;
t366 = t43.*(t364+t365-lambda.*t62.*t81-mu.*t59.*t78).*(1.0./6.0);
t367 = lambda.*t62.*t75;
t368 = mu.*t59.*t77;
t369 = t49.*t50.*(t367+t368).*(1.0./6.0);
t370 = t366+t369;
t372 = t59.*t64.*t81;
t374 = t59.*t64.*t75;
t375 = lambda.*t79.*t81;
t376 = mu.*t79.*t81;
t377 = t375+t376;
t378 = lambda.*t75.*t79;
t379 = mu.*t75.*t79;
t380 = t49.*t50.*(t378+t379).*(1.0./6.0);
t381 = t380-t43.*t377.*(1.0./6.0);
t382 = lambda.*t75.*t78;
t383 = mu.*t75.*t78;
t384 = lambda.*t77.*t81;
t385 = mu.*t77.*t81;
t386 = t382+t383+t384+t385;
t387 = lambda.*(t7-t45-t61+t69+t76-t80).*(t3-t53-t58+t73+t74-t82);
t388 = mu.*(t7-t45-t61+t69+t76-t80).*(t3-t53-t58+t73+t74-t82);
t389 = t387+t388;
t390 = t49.*t50.*t389.*(1.0./6.0);
t391 = t390-t43.*t386.*(1.0./6.0);
t392 = mu.*t77.*t78.*2.0;
t393 = t7-t45-t61+t69+t76-t80;
t394 = mu.*t335;
t395 = t3-t53-t58+t73+t74-t82;
t396 = mu.*t78.*t84;
t397 = mu.*t77.*t84;
t398 = t104+t133;
t399 = t43.*t398.*(1.0./6.0);
t400 = t10.*t49.*t64;
t401 = t105+t134+t400;
t402 = t49.*t50.*t401.*(1.0./6.0);
t403 = t399+t402;
t404 = lambda.*t9.*t49;
t405 = mu.*t10.*t84;
t406 = t404+t405;
t407 = t49.*t50.*t406.*(-1.0./6.0)-lambda.*t26.*t43.*t49.*(1.0./6.0);
t408 = lambda.*t5.*t49;
t409 = mu.*t10.*t83;
t410 = t408+t409;
t411 = t49.*t50.*t410.*(1.0./6.0);
t412 = lambda.*t43.*t49.*t51.*(1.0./6.0);
t413 = t411+t412;
t414 = t198+t249;
t415 = t49.*t64.*t67;
t416 = t199+t250+t415;
t417 = t43.*t414.*(-1.0./6.0)-t49.*t50.*t416.*(1.0./6.0);
t418 = lambda.*t49.*t62;
t419 = mu.*t67.*t84;
t420 = t418+t419;
t421 = t49.*t50.*t420.*(1.0./6.0);
t422 = lambda.*t43.*t49.*t65.*(1.0./6.0);
t423 = t421+t422;
t424 = lambda.*t49.*t59;
t425 = mu.*t67.*t83;
t426 = t424+t425;
t427 = t49.*t50.*t426.*(-1.0./6.0)-lambda.*t43.*t49.*t71.*(1.0./6.0);
t428 = t337+t396;
t429 = t43.*t428.*(1.0./6.0);
t430 = t338+t397-t49.*t64.*t79;
t431 = t429-t49.*t50.*t430.*(1.0./6.0);
t432 = lambda.*t49.*t77;
t433 = t49.*t50.*(t432-mu.*t79.*t84).*(1.0./6.0);
t434 = t433-lambda.*t43.*t49.*t78.*(1.0./6.0);
t435 = lambda.*t49.*t75;
t436 = t435-mu.*t79.*t83;
t437 = lambda.*t43.*t49.*t81.*(1.0./6.0);
t438 = t437-t49.*t50.*t436.*(1.0./6.0);
t439 = lambda.*t10.*t84;
t440 = mu.*t9.*t49;
t441 = t439+t440;
t442 = t49.*t50.*t441.*(-1.0./6.0)-mu.*t26.*t43.*t49.*(1.0./6.0);
t443 = t26.*t64.*t84;
t444 = t104+t443;
t445 = t43.*t444.*(1.0./6.0);
t446 = t9.*t64.*t84;
t447 = t105+t135+t446;
t448 = t49.*t50.*t447.*(1.0./6.0);
t449 = t445+t448;
t450 = mu.*t26.*t83;
t451 = lambda.*t51.*t84;
t452 = t450+t451;
t453 = lambda.*t5.*t84;
t454 = mu.*t9.*t83;
t455 = t453+t454;
t456 = t43.*t452.*(-1.0./6.0)-t49.*t50.*t455.*(1.0./6.0);
t457 = lambda.*t67.*t84;
t458 = mu.*t49.*t62;
t459 = t457+t458;
t460 = t49.*t50.*t459.*(1.0./6.0);
t461 = mu.*t43.*t49.*t65.*(1.0./6.0);
t462 = t460+t461;
t463 = t64.*t65.*t84;
t464 = t198+t463;
t465 = t62.*t64.*t84;
t466 = t199+t251+t465;
t467 = t43.*t464.*(-1.0./6.0)-t49.*t50.*t466.*(1.0./6.0);
t468 = mu.*t65.*t83;
t469 = lambda.*t71.*t84;
t470 = t468+t469;
t471 = t43.*t470.*(1.0./6.0);
t472 = lambda.*t59.*t84;
t473 = mu.*t62.*t83;
t474 = t472+t473;
t475 = t49.*t50.*t474.*(1.0./6.0);
t476 = t471+t475;
t477 = lambda.*t79.*t84;
t478 = t477-mu.*t49.*t77;
t479 = t49.*t50.*t478.*(-1.0./6.0)-mu.*t43.*t49.*t78.*(1.0./6.0);
t480 = t64.*t78.*t84;
t481 = t337+t480;
t482 = t43.*t481.*(1.0./6.0);
t483 = t64.*t77.*t84;
t554 = mu.*t49.*t79;
t484 = t338+t483-t554;
t485 = t482-t49.*t50.*t484.*(1.0./6.0);
t486 = mu.*t78.*t83;
t487 = lambda.*t81.*t84;
t488 = t486+t487;
t489 = lambda.*t75.*t84;
t490 = mu.*t77.*t83;
t491 = t49.*t50.*(t489+t490).*(1.0./6.0);
t492 = t491-t43.*t488.*(1.0./6.0);
t493 = lambda.*t49.*t84;
t494 = mu.*t49.*t84;
t495 = t493+t494;
t496 = t49.*t50.*t495.*(1.0./6.0);
t497 = t83.^2;
t498 = mu.*t497;
t499 = t49.^2;
t500 = t84.^2;
t501 = lambda.*t10.*t83;
t502 = mu.*t5.*t49;
t503 = t501+t502;
t504 = t49.*t50.*t503.*(1.0./6.0);
t505 = mu.*t43.*t49.*t51.*(1.0./6.0);
t506 = t504+t505;
t507 = lambda.*t26.*t83;
t508 = mu.*t51.*t84;
t509 = t507+t508;
t510 = lambda.*t9.*t83;
t511 = mu.*t5.*t84;
t512 = t510+t511;
t513 = t43.*t509.*(-1.0./6.0)-t49.*t50.*t512.*(1.0./6.0);
t514 = t51.*t64.*t83;
t515 = t133+t514;
t516 = t43.*t515.*(1.0./6.0);
t517 = t5.*t64.*t83;
t518 = t134+t135+t517;
t519 = t49.*t50.*t518.*(1.0./6.0);
t520 = t516+t519;
t521 = lambda.*t67.*t83;
t522 = mu.*t49.*t59;
t523 = t521+t522;
t524 = t49.*t50.*t523.*(-1.0./6.0)-mu.*t43.*t49.*t71.*(1.0./6.0);
t525 = lambda.*t65.*t83;
t526 = mu.*t71.*t84;
t527 = t525+t526;
t528 = t43.*t527.*(1.0./6.0);
t529 = lambda.*t62.*t83;
t530 = mu.*t59.*t84;
t531 = t529+t530;
t532 = t49.*t50.*t531.*(1.0./6.0);
t533 = t528+t532;
t534 = t64.*t71.*t83;
t535 = t249+t534;
t536 = t59.*t64.*t83;
t537 = t250+t251+t536;
t538 = t43.*t535.*(-1.0./6.0)-t49.*t50.*t537.*(1.0./6.0);
t539 = lambda.*t79.*t83;
t540 = t539-mu.*t49.*t75;
t541 = t49.*t50.*t540.*(1.0./6.0);
t542 = mu.*t43.*t49.*t81.*(1.0./6.0);
t543 = t541+t542;
t544 = lambda.*t78.*t83;
t545 = mu.*t81.*t84;
t546 = t544+t545;
t547 = lambda.*t77.*t83;
t548 = mu.*t75.*t84;
t549 = t49.*t50.*(t547+t548).*(1.0./6.0);
t550 = t549-t43.*t546.*(1.0./6.0);
t551 = t64.*t81.*t83;
t552 = t396+t551;
t553 = t43.*t552.*(1.0./6.0);
t555 = t64.*t75.*t83;
t556 = lambda.*t49.*t83;
t557 = mu.*t49.*t83;
t558 = t556+t557;
t559 = lambda.*t83.*t84;
t560 = mu.*t83.*t84;
t561 = t559+t560;
t562 = t49.*t50.*t561.*(1.0./6.0);
t563 = mu.*t500;
t564 = mu.*t499;
der_st_dx4 = reshape([t43.*(t94+t124).*(-1.0./6.0)-t49.*t50.*(t96+t125+t64.*t97).*(1.0./6.0),t93,t112,t141,t164,t208,t254,t295,t344,t403,t442,t506,t93,t43.*(t94+t9.*t26.*t64.*2.0).*(-1.0./6.0)-t49.*t50.*(t96+t126+t64.*t98).*(1.0./6.0),t123,t148,t172,t217,t262,t299-t43.*(t102-t296+t297-t26.*t64.*t77).*(1.0./6.0),t352,t407,t449,t513,t112,t123,t43.*(t124+t5.*t51.*t64.*2.0).*(-1.0./6.0)-t49.*t50.*(t125+t126+t64.*t95).*(1.0./6.0),t157,t181,t225,t268,t307,t43.*(t131-t353+t354-t51.*t64.*t75).*(-1.0./6.0)+t49.*t50.*(t132-t355+t356).*(1.0./6.0),t413,t456,t520,t141,t148,t157,t43.*(t191+t244).*(-1.0./6.0)-t49.*t50.*(t193+t245+t64.*t194).*(1.0./6.0),t190,t232,t272,t313,t363,t417,t462,t524,t164,t172,t181,t190,t43.*(t191+t62.*t64.*t65.*2.0).*(-1.0./6.0)-t49.*t50.*(t193+t246+t64.*t195).*(1.0./6.0),t243,t278,t43.*(t196-t314+t315-t64.*t65.*t77).*(1.0./6.0)-t49.*t50.*t317.*(1.0./6.0),t370,t423,t467,t533,t208,t217,t225,t232,t243,t43.*(t244+t59.*t64.*t71.*2.0).*(-1.0./6.0)-t49.*t50.*(t245+t246+t64.*t192).*(1.0./6.0),t286,t324,t43.*(t247-t371+t372-t64.*t71.*t75).*(1.0./6.0)-t49.*t50.*(t248-t373+t374).*(1.0./6.0),t427,t476,t538,t254,t262,t268,t272,t278,t286,t43.*(t333+t392).*(1.0./6.0)-t49.*t50.*(t64.*t335+mu.*t287.^2+mu.*t288.^2).*(1.0./6.0),t332,t381,t431,t479,t543,t295,t299-t43.*(t102+t297-mu.*t51.*t75-t26.*t64.*t77).*(1.0./6.0),t307,t313,t43.*(t196+t315-mu.*t71.*t75-t64.*t65.*t77).*(1.0./6.0)-t49.*t50.*t317.*(1.0./6.0),t324,t332,t43.*(t333+t64.*t77.*t78.*2.0).*(1.0./6.0)-t49.*t50.*(t394+mu.*t334.^2+t64.*t336.^2).*(1.0./6.0),t391,t434,t485,t550,t344,t352,t43.*(t131+t354-mu.*t26.*t77-t51.*t64.*t75).*(-1.0./6.0)+t49.*t50.*(t132+t356-mu.*t10.*t79).*(1.0./6.0),t363,t370,t43.*(t247+t372-mu.*t65.*t77-t64.*t71.*t75).*(1.0./6.0)-t49.*t50.*(t248+t374-mu.*t67.*t79).*(1.0./6.0),t381,t391,t43.*(t392+t64.*t75.*t81.*2.0).*(1.0./6.0)-t49.*t50.*(t394+mu.*t393.^2+t64.*t395.^2).*(1.0./6.0),t438,t492,t553-t49.*t50.*(t397-t554+t555).*(1.0./6.0),t403,t407,t413,t417,t423,t427,t431,t434,t438,t49.*t50.*(t498+t563+t64.*t499).*(-1.0./6.0),t496,t49.*t50.*t558.*(-1.0./6.0),t442,t449,t456,t462,t467,t476,t479,t485,t492,t496,t49.*t50.*(t498+t564+t64.*t500).*(-1.0./6.0),t562,t506,t513,t520,t524,t533,t538,t543,t550,t553-t49.*t50.*(t397+t555-mu.*t49.*t79).*(1.0./6.0),t49.*t50.*t558.*(-1.0./6.0),t562,t49.*t50.*(t563+t564+t64.*t497).*(-1.0./6.0)],[12, 12]);
