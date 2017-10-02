function der_st_dz3 = der_stima_dz3(lambda,mu,V)
%DER_STIMA_DZ3
%    DER_ST_DZ3 = DER_STIMA_DZ3(LAMBDA,MU,X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    02-Jun-2015 16:44:13
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
t13 = x2.*y4;
t14 = x4.*y2;
t53 = x2.*y3;
t54 = x3.*y2;
t55 = x3.*y4;
t56 = x4.*y3;
t2 = -t13+t14+t53-t54+t55-t56;
t3 = x2.*z3;
t4 = x4.*z2;
t5 = x3.*z4;
t29 = x3.*z2;
t30 = x2.*z4;
t31 = x4.*z3;
t6 = t3+t4+t5-t29-t30-t31;
t7 = mu.*2.0;
t8 = lambda+t7;
t9 = y2.*z3;
t10 = y4.*z2;
t11 = y3.*z4;
t33 = y3.*z2;
t34 = y2.*z4;
t35 = y4.*z3;
t12 = t9+t10+t11-t33-t34-t35;
t15 = x1.*y2.*z3;
t16 = x2.*y3.*z1;
t17 = x3.*y1.*z2;
t18 = x1.*y4.*z2;
t19 = x2.*y1.*z4;
t20 = x4.*y2.*z1;
t21 = x1.*y3.*z4;
t22 = x3.*y4.*z1;
t23 = x4.*y1.*z3;
t24 = x2.*y4.*z3;
t25 = x3.*y2.*z4;
t26 = x4.*y3.*z2;
t36 = x1.*y3.*z2;
t37 = x2.*y1.*z3;
t38 = x3.*y2.*z1;
t39 = x1.*y2.*z4;
t40 = x2.*y4.*z1;
t41 = x4.*y1.*z2;
t42 = x1.*y4.*z3;
t43 = x3.*y1.*z4;
t44 = x4.*y3.*z1;
t45 = x2.*y3.*z4;
t46 = x3.*y4.*z2;
t47 = x4.*y2.*z3;
t27 = t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25+t26-t36-t37-t38-t39-t40-t41-t42-t43-t44-t45-t46-t47;
t28 = y2-y4;
t32 = x2-x4;
t48 = 1.0./t27;
t49 = x1.*y2;
t50 = x4.*y1;
t57 = x2.*y1;
t58 = x1.*y4;
t51 = t13-t14+t49+t50-t57-t58;
t52 = 1.0./t27.^2;
t59 = x1.*z3;
t60 = x4.*z1;
t65 = x3.*z1;
t66 = x1.*z4;
t61 = t5-t31+t59+t60-t65-t66;
t62 = y1.*z3;
t63 = y4.*z1;
t69 = y3.*z1;
t70 = y1.*z4;
t64 = t11-t35+t62+t63-t69-t70;
t67 = y1-y4;
t68 = x1-x4;
t71 = x1.*y3;
t73 = x3.*y1;
t72 = t50+t55-t56-t58+t71-t73;
t74 = x2.*z1;
t78 = x1.*z2;
t75 = t4-t30-t60+t66+t74-t78;
t76 = y2.*z1;
t79 = y1.*z2;
t77 = t10-t34-t63+t70+t76-t79;
t80 = t3-t29-t59+t65-t74+t78;
t81 = t9-t33-t62+t69-t76+t79;
t82 = y1-y2;
t83 = x1-x2;
t84 = t49+t53-t54-t57-t71+t73;
t85 = lambda.*t6.*t28;
t86 = mu.*t6.*t28;
t87 = lambda.*t12.*t32;
t88 = mu.*t12.*t32;
t89 = t85+t86+t87+t88;
t90 = t48.*t89.*(1.0./6.0);
t91 = lambda.*t6.*t12;
t92 = mu.*t6.*t12;
t93 = t91+t92;
t94 = t90-t51.*t52.*t93.*(1.0./6.0);
t95 = t2.^2;
t96 = mu.*t95;
t97 = t12.^2;
t98 = t6.^2;
t99 = mu.*t2.*t72;
t100 = mu.*t2.*t84;
t101 = lambda.*t2.*t28;
t102 = mu.*t2.*t28;
t103 = t101+t102;
t104 = lambda.*t2.*t12;
t105 = mu.*t2.*t12;
t106 = t104+t105;
t107 = t51.*t52.*t106.*(1.0./6.0);
t108 = t107-t48.*t103.*(1.0./6.0);
t109 = lambda.*t2.*t32;
t110 = mu.*t2.*t32;
t111 = t109+t110;
t112 = t48.*t111.*(1.0./6.0);
t113 = lambda.*t2.*t6;
t114 = mu.*t2.*t6;
t115 = t113+t114;
t116 = t112-t51.*t52.*t115.*(1.0./6.0);
t117 = mu.*t6.*t32.*2.0;
t118 = mu.*t12.*t28.*2.0;
t119 = mu.*t98;
t120 = mu.*t97;
t121 = mu.*t32.*t61;
t122 = mu.*t6.*t68;
t123 = mu.*t28.*t64;
t124 = mu.*t12.*t67;
t125 = mu.*t6.*t61;
t126 = mu.*t12.*t64;
t127 = mu.*t32.*t75;
t128 = mu.*t28.*t77;
t129 = mu.*t6.*t75;
t130 = mu.*t12.*t77;
t131 = mu.*t6.*t83;
t132 = mu.*t12.*t82;
t133 = mu.*t6.*t80;
t134 = mu.*t12.*t81;
t135 = t8.*t28.*t64;
t136 = t8.*t12.*t67;
t137 = t121+t122+t135+t136;
t138 = t48.*t137.*(1.0./6.0);
t139 = t8.*t12.*t64;
t140 = t99+t125+t139;
t141 = t138-t51.*t52.*t140.*(1.0./6.0);
t142 = lambda.*t6.*t67;
t143 = mu.*t28.*t61;
t144 = lambda.*t32.*t64;
t145 = mu.*t12.*t68;
t146 = t142+t143+t144+t145;
t147 = lambda.*t6.*t64;
t148 = mu.*t12.*t61;
t149 = t147+t148;
t150 = t51.*t52.*t149.*(1.0./6.0);
t151 = t150-t48.*t146.*(1.0./6.0);
t152 = lambda.*t2.*t67;
t153 = mu.*t28.*t72;
t154 = t152+t153;
t155 = t48.*t154.*(1.0./6.0);
t156 = lambda.*t2.*t64;
t157 = mu.*t12.*t72;
t158 = t156+t157;
t159 = t155-t51.*t52.*t158.*(1.0./6.0);
t160 = lambda.*t28.*t61;
t161 = mu.*t6.*t67;
t162 = lambda.*t12.*t68;
t163 = mu.*t32.*t64;
t164 = t160+t161+t162+t163;
t165 = lambda.*t12.*t61;
t166 = mu.*t6.*t64;
t167 = t165+t166;
t168 = t51.*t52.*t167.*(1.0./6.0);
t169 = t168-t48.*t164.*(1.0./6.0);
t170 = t8.*t32.*t61;
t171 = t6.*t8.*t68;
t172 = t123+t124+t170+t171;
t173 = t48.*t172.*(1.0./6.0);
t174 = t6.*t8.*t61;
t175 = t99+t126+t174;
t176 = t173-t51.*t52.*t175.*(1.0./6.0);
t177 = lambda.*t2.*t68;
t178 = mu.*t32.*t72;
t179 = t177+t178;
t180 = lambda.*t2.*t61;
t181 = mu.*t6.*t72;
t182 = t180+t181;
t183 = t51.*t52.*t182.*(1.0./6.0);
t184 = t183-t48.*t179.*(1.0./6.0);
t185 = lambda.*t61.*t67;
t186 = mu.*t61.*t67;
t187 = lambda.*t64.*t68;
t188 = mu.*t64.*t68;
t189 = t185+t186+t187+t188;
t190 = t48.*t189.*(1.0./6.0);
t191 = lambda.*t61.*t64;
t192 = mu.*t61.*t64;
t193 = t191+t192;
t194 = t190-t51.*t52.*t193.*(1.0./6.0);
t195 = t72.^2;
t196 = mu.*t195;
t197 = t64.^2;
t198 = t61.^2;
t199 = mu.*t72.*t84;
t200 = lambda.*t28.*t72;
t201 = mu.*t2.*t67;
t202 = t200+t201;
t203 = t48.*t202.*(1.0./6.0);
t204 = lambda.*t12.*t72;
t205 = mu.*t2.*t64;
t206 = t204+t205;
t207 = t203-t51.*t52.*t206.*(1.0./6.0);
t208 = lambda.*t32.*t72;
t209 = mu.*t2.*t68;
t210 = t208+t209;
t211 = lambda.*t6.*t72;
t212 = mu.*t2.*t61;
t213 = t211+t212;
t214 = t51.*t52.*t213.*(1.0./6.0);
t215 = t214-t48.*t210.*(1.0./6.0);
t216 = t121+t122+t123+t124;
t217 = t48.*t216.*(1.0./6.0);
t218 = t2.*t8.*t72;
t219 = t125+t126+t218;
t220 = t217-t51.*t52.*t219.*(1.0./6.0);
t221 = lambda.*t67.*t72;
t222 = mu.*t67.*t72;
t223 = t221+t222;
t224 = lambda.*t64.*t72;
t225 = mu.*t64.*t72;
t226 = t224+t225;
t227 = t51.*t52.*t226.*(1.0./6.0);
t228 = t227-t48.*t223.*(1.0./6.0);
t229 = lambda.*t68.*t72;
t230 = mu.*t68.*t72;
t231 = t229+t230;
t232 = t48.*t231.*(1.0./6.0);
t233 = lambda.*t61.*t72;
t234 = mu.*t61.*t72;
t235 = t233+t234;
t236 = t232-t51.*t52.*t235.*(1.0./6.0);
t237 = mu.*t61.*t68.*2.0;
t238 = mu.*t64.*t67.*2.0;
t239 = mu.*t198;
t240 = mu.*t197;
t241 = mu.*t68.*t75;
t242 = mu.*t67.*t77;
t243 = mu.*t61.*t75;
t244 = mu.*t64.*t77;
t245 = mu.*t61.*t83;
t246 = mu.*t64.*t82;
t247 = mu.*t61.*t80;
t248 = mu.*t64.*t81;
t249 = t8.*t28.*t77;
t250 = t48.*(t127+t249).*(1.0./6.0);
t251 = t8.*t12.*t77;
t291 = mu.*t2.*t51;
t252 = t129+t251-t291;
t253 = t250-t51.*t52.*t252.*(1.0./6.0);
t254 = mu.*t28.*t75;
t255 = lambda.*t32.*t77;
t256 = t254+t255;
t257 = lambda.*t6.*t77;
t258 = mu.*t12.*t75;
t259 = t51.*t52.*(t257+t258).*(1.0./6.0);
t260 = t259-t48.*t256.*(1.0./6.0);
t261 = lambda.*t2.*t77;
t262 = t261-mu.*t12.*t51;
t263 = t51.*t52.*t262.*(-1.0./6.0)-mu.*t28.*t48.*t51.*(1.0./6.0);
t264 = t8.*t67.*t77;
t265 = t241+t264;
t266 = t8.*t64.*t77;
t306 = mu.*t51.*t72;
t267 = t51.*t52.*(t243+t266-t306).*(1.0./6.0);
t268 = t267-t48.*t265.*(1.0./6.0);
t269 = mu.*t67.*t75;
t270 = lambda.*t68.*t77;
t271 = t48.*(t269+t270).*(1.0./6.0);
t272 = lambda.*t61.*t77;
t273 = mu.*t64.*t75;
t274 = t272+t273;
t275 = t271-t51.*t52.*t274.*(1.0./6.0);
t276 = lambda.*t72.*t77;
t277 = t51.*t52.*(t276-mu.*t51.*t64).*(1.0./6.0);
t278 = mu.*t48.*t51.*t67.*(1.0./6.0);
t279 = t277+t278;
t280 = t4-t30-t60+t66+t74-t78;
t281 = t10-t34-t63+t70+t76-t79;
t282 = lambda.*t28.*t75;
t283 = mu.*t32.*t77;
t284 = t282+t283;
t285 = lambda.*t12.*t75;
t286 = mu.*t6.*t77;
t287 = t51.*t52.*(t285+t286).*(1.0./6.0);
t288 = t287-t48.*t284.*(1.0./6.0);
t289 = t8.*t32.*t75;
t290 = t48.*(t128+t289).*(1.0./6.0);
t292 = t6.*t8.*t75;
t293 = lambda.*t2.*t75;
t294 = t51.*t52.*(t293-mu.*t6.*t51).*(1.0./6.0);
t295 = mu.*t32.*t48.*t51.*(1.0./6.0);
t296 = t294+t295;
t297 = lambda.*t67.*t75;
t298 = mu.*t68.*t77;
t299 = t48.*(t297+t298).*(1.0./6.0);
t300 = lambda.*t64.*t75;
t301 = mu.*t61.*t77;
t302 = t300+t301;
t303 = t299-t51.*t52.*t302.*(1.0./6.0);
t304 = t8.*t68.*t75;
t305 = t242+t304;
t307 = t8.*t61.*t75;
t308 = lambda.*t72.*t75;
t309 = t308-mu.*t51.*t61;
t310 = t51.*t52.*t309.*(-1.0./6.0)-mu.*t48.*t51.*t68.*(1.0./6.0);
t311 = lambda.*t75.*t77;
t312 = mu.*t75.*t77;
t313 = t311+t312;
t314 = t51.^2;
t315 = mu.*t314;
t316 = t10-t34-t63+t70+t76-t79;
t317 = t4-t30-t60+t66+t74-t78;
t318 = lambda.*t12.*t51;
t319 = t318-mu.*t2.*t77;
t320 = t51.*t52.*t319.*(1.0./6.0);
t321 = t320-lambda.*t28.*t48.*t51.*(1.0./6.0);
t322 = lambda.*t6.*t51;
t323 = t322-mu.*t2.*t75;
t324 = lambda.*t32.*t48.*t51.*(1.0./6.0);
t325 = t324-t51.*t52.*t323.*(1.0./6.0);
t326 = t48.*(t127+t128).*(1.0./6.0);
t327 = t129+t130-t2.*t8.*t51;
t328 = t326-t51.*t52.*t327.*(1.0./6.0);
t329 = lambda.*t51.*t64;
t330 = t329-mu.*t72.*t77;
t331 = lambda.*t48.*t51.*t67.*(1.0./6.0);
t332 = t331-t51.*t52.*t330.*(1.0./6.0);
t333 = lambda.*t51.*t61;
t334 = t333-mu.*t72.*t75;
t335 = t51.*t52.*t334.*(1.0./6.0);
t336 = t335-lambda.*t48.*t51.*t68.*(1.0./6.0);
t337 = t241+t242;
t338 = t51.*t52.*(t243+t244-t8.*t51.*t72).*(1.0./6.0);
t339 = t338-t48.*t337.*(1.0./6.0);
t340 = lambda.*t51.*t77;
t341 = mu.*t51.*t77;
t342 = t340+t341;
t343 = lambda.*t51.*t75;
t344 = mu.*t51.*t75;
t345 = t51.*t52.*(t343+t344).*(1.0./6.0);
t346 = t4-t30-t60+t66+t74-t78;
t347 = t10-t34-t63+t70+t76-t79;
t348 = mu.*t75.*t83;
t349 = mu.*t77.*t82;
t350 = mu.*t75.*t80;
t351 = mu.*t77.*t81;
t352 = t8.*t12.*t82;
t491 = mu.*t32.*t80;
t353 = t131+t352-t491-t8.*t28.*t81;
t354 = t8.*t12.*t81;
t355 = t100+t133+t354;
t356 = t48.*t353.*(-1.0./6.0)-t51.*t52.*t355.*(1.0./6.0);
t357 = lambda.*t6.*t82;
t358 = mu.*t12.*t83;
t359 = t357+t358-lambda.*t32.*t81-mu.*t28.*t80;
t360 = t48.*t359.*(1.0./6.0);
t361 = lambda.*t6.*t81;
t362 = mu.*t12.*t80;
t363 = t361+t362;
t364 = t51.*t52.*t363.*(1.0./6.0);
t365 = t360+t364;
t366 = lambda.*t2.*t82;
t367 = t366-mu.*t28.*t84;
t368 = lambda.*t2.*t81;
t369 = mu.*t12.*t84;
t370 = t368+t369;
t371 = t48.*t367.*(-1.0./6.0)-t51.*t52.*t370.*(1.0./6.0);
t372 = t8.*t64.*t82;
t509 = mu.*t68.*t80;
t373 = t48.*(t245+t372-t509-t8.*t67.*t81).*(1.0./6.0);
t374 = t8.*t64.*t81;
t375 = t199+t247+t374;
t376 = t51.*t52.*t375.*(1.0./6.0);
t377 = t373+t376;
t378 = lambda.*t61.*t82;
t379 = mu.*t64.*t83;
t380 = t378+t379-lambda.*t68.*t81-mu.*t67.*t80;
t381 = lambda.*t61.*t81;
t382 = mu.*t64.*t80;
t383 = t381+t382;
t384 = t48.*t380.*(-1.0./6.0)-t51.*t52.*t383.*(1.0./6.0);
t385 = lambda.*t72.*t82;
t386 = t385-mu.*t67.*t84;
t387 = t48.*t386.*(1.0./6.0);
t388 = lambda.*t72.*t81;
t389 = mu.*t64.*t84;
t390 = t388+t389;
t391 = t51.*t52.*t390.*(1.0./6.0);
t392 = t387+t391;
t393 = t8.*t77.*t82;
t394 = t48.*(t348+t393).*(1.0./6.0);
t395 = t8.*t77.*t81;
t457 = mu.*t51.*t84;
t396 = t51.*t52.*(t350+t395-t457).*(1.0./6.0);
t397 = t394+t396;
t398 = lambda.*t75.*t82;
t399 = mu.*t77.*t83;
t400 = t398+t399;
t401 = lambda.*t75.*t81;
t402 = mu.*t77.*t80;
t403 = t401+t402;
t404 = t48.*t400.*(-1.0./6.0)-t51.*t52.*t403.*(1.0./6.0);
t405 = lambda.*t51.*t81;
t406 = t405-mu.*t77.*t84;
t407 = t51.*t52.*t406.*(-1.0./6.0)-lambda.*t48.*t51.*t82.*(1.0./6.0);
t408 = mu.*t6.*t82;
t409 = lambda.*t12.*t83;
t410 = t48.*(t408+t409-lambda.*t28.*t80-mu.*t32.*t81).*(1.0./6.0);
t411 = lambda.*t12.*t80;
t412 = mu.*t6.*t81;
t413 = t411+t412;
t414 = t51.*t52.*t413.*(1.0./6.0);
t415 = t410+t414;
t416 = t6.*t8.*t83;
t492 = mu.*t28.*t81;
t417 = t132+t416-t492-t8.*t32.*t80;
t418 = t6.*t8.*t80;
t419 = t100+t134+t418;
t420 = t48.*t417.*(-1.0./6.0)-t51.*t52.*t419.*(1.0./6.0);
t421 = lambda.*t2.*t83;
t422 = t421-mu.*t32.*t84;
t423 = t48.*t422.*(1.0./6.0);
t424 = lambda.*t2.*t80;
t425 = mu.*t6.*t84;
t426 = t424+t425;
t427 = t51.*t52.*t426.*(1.0./6.0);
t428 = t423+t427;
t429 = mu.*t61.*t82;
t430 = lambda.*t64.*t83;
t431 = t429+t430-lambda.*t67.*t80-mu.*t68.*t81;
t432 = lambda.*t64.*t80;
t433 = mu.*t61.*t81;
t434 = t432+t433;
t435 = t48.*t431.*(-1.0./6.0)-t51.*t52.*t434.*(1.0./6.0);
t436 = t8.*t61.*t83;
t510 = mu.*t67.*t81;
t437 = t48.*(t246+t436-t510-t8.*t68.*t80).*(1.0./6.0);
t438 = t8.*t61.*t80;
t439 = t199+t248+t438;
t440 = t51.*t52.*t439.*(1.0./6.0);
t441 = t437+t440;
t442 = lambda.*t72.*t83;
t443 = t442-mu.*t68.*t84;
t444 = lambda.*t72.*t80;
t445 = mu.*t61.*t84;
t446 = t444+t445;
t447 = t48.*t443.*(-1.0./6.0)-t51.*t52.*t446.*(1.0./6.0);
t448 = mu.*t75.*t82;
t449 = lambda.*t77.*t83;
t450 = t448+t449;
t451 = lambda.*t77.*t80;
t452 = mu.*t75.*t81;
t453 = t451+t452;
t454 = t48.*t450.*(-1.0./6.0)-t51.*t52.*t453.*(1.0./6.0);
t455 = t8.*t75.*t83;
t456 = t48.*(t349+t455).*(1.0./6.0);
t458 = t8.*t75.*t80;
t459 = lambda.*t51.*t80;
t460 = t459-mu.*t75.*t84;
t461 = t51.*t52.*t460.*(1.0./6.0);
t462 = lambda.*t48.*t51.*t83.*(1.0./6.0);
t463 = t461+t462;
t464 = lambda.*t80.*t82;
t465 = mu.*t80.*t82;
t466 = lambda.*t81.*t83;
t467 = mu.*t81.*t83;
t468 = t464+t465+t466+t467;
t469 = lambda.*t80.*t81;
t470 = mu.*t80.*t81;
t471 = t469+t470;
t472 = t48.*t468.*(-1.0./6.0)-t51.*t52.*t471.*(1.0./6.0);
t473 = t84.^2;
t474 = mu.*t473;
t475 = t81.^2;
t476 = t80.^2;
t477 = lambda.*t28.*t84;
t478 = t477-mu.*t2.*t82;
t479 = t48.*t478.*(1.0./6.0);
t480 = lambda.*t12.*t84;
t481 = mu.*t2.*t81;
t482 = t480+t481;
t483 = t479-t51.*t52.*t482.*(1.0./6.0);
t484 = lambda.*t32.*t84;
t485 = t484-mu.*t2.*t83;
t486 = lambda.*t6.*t84;
t487 = mu.*t2.*t80;
t488 = t486+t487;
t489 = t51.*t52.*t488.*(1.0./6.0);
t490 = t489-t48.*t485.*(1.0./6.0);
t493 = t2.*t8.*t84;
t494 = t133+t134+t493;
t495 = lambda.*t67.*t84;
t496 = t495-mu.*t72.*t82;
t497 = lambda.*t64.*t84;
t498 = mu.*t72.*t81;
t499 = t497+t498;
t500 = t51.*t52.*t499.*(1.0./6.0);
t501 = t500-t48.*t496.*(1.0./6.0);
t502 = lambda.*t68.*t84;
t503 = t502-mu.*t72.*t83;
t504 = t48.*t503.*(1.0./6.0);
t505 = lambda.*t61.*t84;
t506 = mu.*t72.*t80;
t507 = t505+t506;
t508 = t504-t51.*t52.*t507.*(1.0./6.0);
t511 = t8.*t72.*t84;
t512 = t247+t248+t511;
t513 = t51.*t52.*t512.*(1.0./6.0);
t514 = lambda.*t77.*t84;
t515 = t51.*t52.*(t514-mu.*t51.*t81).*(1.0./6.0);
t516 = t515-mu.*t48.*t51.*t82.*(1.0./6.0);
t517 = lambda.*t75.*t84;
t518 = t517-mu.*t51.*t80;
t519 = mu.*t48.*t51.*t83.*(1.0./6.0);
t520 = t519-t51.*t52.*t518.*(1.0./6.0);
t521 = t48.*(t348+t349).*(1.0./6.0);
t522 = t51.*t52.*(t350+t351-t8.*t51.*t84).*(1.0./6.0);
t523 = t521+t522;
t524 = lambda.*t82.*t84;
t525 = mu.*t82.*t84;
t526 = t524+t525;
t527 = t48.*t526.*(1.0./6.0);
t528 = lambda.*t81.*t84;
t529 = mu.*t81.*t84;
t530 = t528+t529;
t531 = t51.*t52.*t530.*(1.0./6.0);
t532 = t527+t531;
t533 = lambda.*t83.*t84;
t534 = mu.*t83.*t84;
t535 = t533+t534;
t536 = lambda.*t80.*t84;
t537 = mu.*t80.*t84;
t538 = t536+t537;
t539 = t48.*t535.*(-1.0./6.0)-t51.*t52.*t538.*(1.0./6.0);
t540 = mu.*t80.*t83.*2.0;
t541 = mu.*t81.*t82.*2.0;
t542 = mu.*t476;
t543 = mu.*t475;
der_st_dz3 = reshape([t48.*(t117+t8.*t12.*t28.*2.0).*(-1.0./6.0)+t51.*t52.*(t96+t119+t8.*t97).*(1.0./6.0),t94,t108,t141,t169,t207,t253,t288,t321,t356,t415,t483,t94,t48.*(t118+t6.*t8.*t32.*2.0).*(-1.0./6.0)+t51.*t52.*(t96+t120+t8.*t98).*(1.0./6.0),t116,t151,t176,t215,t260,t290-t51.*t52.*(t130-t291+t292).*(1.0./6.0),t325,t365,t420,t490,t108,t116,t48.*(t117+t118).*(-1.0./6.0)+t51.*t52.*(t119+t120+t8.*t95).*(1.0./6.0),t159,t184,t220,t263,t296,t328,t371,t428,t48.*(t131+t132-t491-t492).*(-1.0./6.0)-t51.*t52.*t494.*(1.0./6.0),t141,t151,t159,t48.*(t237+t8.*t64.*t67.*2.0).*(-1.0./6.0)+t51.*t52.*(t196+t239+t8.*t197).*(1.0./6.0),t194,t228,t268,t303,t332,t377,t435,t501,t169,t176,t184,t194,t48.*(t238+t8.*t61.*t68.*2.0).*(-1.0./6.0)+t51.*t52.*(t196+t240+t8.*t198).*(1.0./6.0),t236,t275,t48.*t305.*(-1.0./6.0)+t51.*t52.*(t244-t306+t307).*(1.0./6.0),t336,t384,t441,t508,t207,t215,t220,t228,t236,t48.*(t237+t238).*(-1.0./6.0)+t51.*t52.*(t239+t240+t8.*t195).*(1.0./6.0),t279,t310,t339,t392,t447,t513+t48.*(t245+t246-t509-t510).*(1.0./6.0),t253,t260,t263,t268,t275,t279,t51.*t52.*(t315+mu.*t280.^2+t8.*t281.^2).*(1.0./6.0),t51.*t52.*t313.*(-1.0./6.0),t51.*t52.*t342.*(-1.0./6.0),t397,t454,t516,t288,t290-t51.*t52.*(t130+t292-mu.*t2.*t51).*(1.0./6.0),t296,t303,t48.*t305.*(-1.0./6.0)+t51.*t52.*(t244+t307-mu.*t51.*t72).*(1.0./6.0),t310,t51.*t52.*t313.*(-1.0./6.0),t51.*t52.*(t315+mu.*t316.^2+t8.*t317.^2).*(1.0./6.0),t345,t404,t456+t51.*t52.*(t351-t457+t458).*(1.0./6.0),t520,t321,t325,t328,t332,t336,t339,t51.*t52.*t342.*(-1.0./6.0),t345,t51.*t52.*(t8.*t314+mu.*t346.^2+mu.*t347.^2).*(1.0./6.0),t407,t463,t523,t356,t365,t371,t377,t384,t392,t397,t404,t407,t48.*(t540+t8.*t81.*t82.*2.0).*(1.0./6.0)+t51.*t52.*(t474+t542+t8.*t475).*(1.0./6.0),t472,t532,t415,t420,t428,t435,t441,t447,t454,t456+t51.*t52.*(t351+t458-mu.*t51.*t84).*(1.0./6.0),t463,t472,t48.*(t541+t8.*t80.*t83.*2.0).*(1.0./6.0)+t51.*t52.*(t474+t543+t8.*t476).*(1.0./6.0),t539,t483,t490,t48.*(t131+t132-mu.*t28.*t81-mu.*t32.*t80).*(-1.0./6.0)-t51.*t52.*t494.*(1.0./6.0),t501,t508,t513+t48.*(t245+t246-mu.*t67.*t81-mu.*t68.*t80).*(1.0./6.0),t516,t520,t523,t532,t539,t48.*(t540+t541).*(1.0./6.0)+t51.*t52.*(t542+t543+t8.*t473).*(1.0./6.0)],[12, 12]);