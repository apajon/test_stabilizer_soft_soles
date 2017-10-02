function der_st_dz4 = der_stima_dz4(lambda,mu,V)
%DER_STIMA_DZ4
%    DER_ST_DZ4 = DER_STIMA_DZ4(LAMBDA,MU,X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    02-Jun-2015 16:46:56
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
t13 = x2.*y3;
t14 = x3.*y2;
t53 = x2.*y4;
t54 = x4.*y2;
t55 = x3.*y4;
t56 = x4.*y3;
t2 = t13-t14-t53+t54+t55-t56;
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
t28 = y2-y3;
t32 = x2-x3;
t48 = 1.0./t27;
t49 = x1.*y2;
t50 = x3.*y1;
t57 = x2.*y1;
t58 = x1.*y3;
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
t67 = y1-y3;
t68 = x1-x3;
t71 = x1.*y4;
t73 = x4.*y1;
t72 = t50-t55+t56-t58+t71-t73;
t74 = x2.*z1;
t78 = x1.*z2;
t75 = t4-t30-t60+t66+t74-t78;
t76 = y2.*z1;
t81 = y1.*z2;
t77 = t10-t34-t63+t70+t76-t81;
t79 = y1-y2;
t80 = x1-x2;
t82 = t49+t53-t54-t57-t71+t73;
t83 = t3-t29-t59+t65-t74+t78;
t84 = t9-t33-t62+t69-t76+t81;
t85 = lambda.*t6.*t28;
t86 = mu.*t6.*t28;
t87 = lambda.*t12.*t32;
t88 = mu.*t12.*t32;
t89 = t85+t86+t87+t88;
t90 = lambda.*t6.*t12;
t91 = mu.*t6.*t12;
t92 = t90+t91;
t93 = t51.*t52.*t92.*(1.0./6.0);
t94 = t93-t48.*t89.*(1.0./6.0);
t95 = t2.^2;
t96 = mu.*t95;
t97 = t12.^2;
t98 = t6.^2;
t99 = mu.*t2.*t51;
t100 = lambda.*t2.*t28;
t101 = mu.*t2.*t28;
t102 = t100+t101;
t103 = t48.*t102.*(1.0./6.0);
t104 = lambda.*t2.*t12;
t105 = mu.*t2.*t12;
t106 = t104+t105;
t107 = t103-t51.*t52.*t106.*(1.0./6.0);
t108 = lambda.*t2.*t32;
t109 = mu.*t2.*t32;
t110 = t108+t109;
t111 = lambda.*t2.*t6;
t112 = mu.*t2.*t6;
t113 = t111+t112;
t114 = t51.*t52.*t113.*(1.0./6.0);
t115 = t114-t48.*t110.*(1.0./6.0);
t116 = mu.*t6.*t32.*2.0;
t117 = mu.*t12.*t28.*2.0;
t118 = mu.*t98;
t119 = mu.*t97;
t120 = mu.*t32.*t61;
t121 = mu.*t6.*t68;
t122 = mu.*t28.*t64;
t123 = mu.*t12.*t67;
t124 = mu.*t6.*t61;
t125 = mu.*t12.*t64;
t126 = mu.*t6.*t80;
t127 = mu.*t12.*t79;
t128 = mu.*t6.*t75;
t129 = mu.*t12.*t77;
t130 = mu.*t32.*t83;
t131 = mu.*t28.*t84;
t132 = mu.*t6.*t83;
t133 = mu.*t12.*t84;
t134 = t8.*t28.*t64;
t135 = t8.*t12.*t67;
t136 = t120+t121+t134+t135;
t137 = t8.*t12.*t64;
t171 = mu.*t2.*t72;
t138 = t124+t137-t171;
t139 = t51.*t52.*t138.*(1.0./6.0);
t140 = t139-t48.*t136.*(1.0./6.0);
t141 = lambda.*t6.*t67;
t142 = mu.*t28.*t61;
t143 = lambda.*t32.*t64;
t144 = mu.*t12.*t68;
t145 = t141+t142+t143+t144;
t146 = t48.*t145.*(1.0./6.0);
t147 = lambda.*t6.*t64;
t148 = mu.*t12.*t61;
t149 = t147+t148;
t150 = t146-t51.*t52.*t149.*(1.0./6.0);
t151 = lambda.*t2.*t67;
t152 = t151-mu.*t28.*t72;
t153 = lambda.*t2.*t64;
t154 = t153-mu.*t12.*t72;
t155 = t51.*t52.*t154.*(1.0./6.0);
t156 = t155-t48.*t152.*(1.0./6.0);
t157 = t50-t55+t56-t58+t71-t73;
t158 = lambda.*t28.*t61;
t159 = mu.*t6.*t67;
t160 = lambda.*t12.*t68;
t161 = mu.*t32.*t64;
t162 = t158+t159+t160+t161;
t163 = t48.*t162.*(1.0./6.0);
t164 = lambda.*t12.*t61;
t165 = mu.*t6.*t64;
t166 = t164+t165;
t167 = t163-t51.*t52.*t166.*(1.0./6.0);
t168 = t8.*t32.*t61;
t169 = t6.*t8.*t68;
t170 = t122+t123+t168+t169;
t172 = t6.*t8.*t61;
t173 = lambda.*t2.*t68;
t174 = t173-mu.*t32.*t72;
t175 = t48.*t174.*(1.0./6.0);
t176 = lambda.*t2.*t61;
t177 = t176-mu.*t6.*t72;
t178 = t175-t51.*t52.*t177.*(1.0./6.0);
t179 = lambda.*t61.*t67;
t180 = mu.*t61.*t67;
t181 = lambda.*t64.*t68;
t182 = mu.*t64.*t68;
t183 = t179+t180+t181+t182;
t184 = lambda.*t61.*t64;
t185 = mu.*t61.*t64;
t186 = t184+t185;
t187 = t51.*t52.*t186.*(1.0./6.0);
t188 = t187-t48.*t183.*(1.0./6.0);
t189 = t50-t55+t56-t58+t71-t73;
t190 = t64.^2;
t191 = t61.^2;
t192 = mu.*t72.*t82;
t193 = lambda.*t28.*t72;
t194 = t48.*(t193-mu.*t2.*t67).*(1.0./6.0);
t195 = lambda.*t12.*t72;
t196 = t195-mu.*t2.*t64;
t197 = t194-t51.*t52.*t196.*(1.0./6.0);
t198 = lambda.*t32.*t72;
t199 = t198-mu.*t2.*t68;
t200 = lambda.*t6.*t72;
t201 = t51.*t52.*(t200-mu.*t2.*t61).*(1.0./6.0);
t202 = t201-t48.*t199.*(1.0./6.0);
t203 = t120+t121+t122+t123;
t204 = t124+t125-t2.*t8.*t72;
t205 = t51.*t52.*t204.*(1.0./6.0);
t206 = t205-t48.*t203.*(1.0./6.0);
t207 = lambda.*t67.*t72;
t208 = mu.*t67.*t72;
t209 = t207+t208;
t210 = lambda.*t64.*t72;
t211 = mu.*t64.*t72;
t212 = t51.*t52.*(t210+t211).*(1.0./6.0);
t213 = t212-t48.*t209.*(1.0./6.0);
t214 = lambda.*t68.*t72;
t215 = mu.*t68.*t72;
t216 = t48.*(t214+t215).*(1.0./6.0);
t217 = lambda.*t61.*t72;
t218 = mu.*t61.*t72;
t219 = t217+t218;
t220 = t216-t51.*t52.*t219.*(1.0./6.0);
t221 = mu.*t61.*t68.*2.0;
t222 = mu.*t64.*t67.*2.0;
t223 = mu.*t191;
t224 = mu.*t190;
t225 = t50-t55+t56-t58+t71-t73;
t226 = mu.*t61.*t80;
t227 = mu.*t64.*t79;
t228 = mu.*t61.*t75;
t229 = mu.*t64.*t77;
t230 = mu.*t68.*t83;
t231 = mu.*t67.*t84;
t232 = mu.*t61.*t83;
t233 = mu.*t64.*t84;
t234 = t8.*t12.*t79;
t340 = mu.*t32.*t75;
t235 = t126+t234-t340-t8.*t28.*t77;
t236 = t48.*t235.*(1.0./6.0);
t237 = t8.*t12.*t77;
t286 = mu.*t2.*t82;
t238 = t51.*t52.*(t128+t237-t286).*(1.0./6.0);
t239 = t236+t238;
t240 = lambda.*t6.*t79;
t241 = mu.*t12.*t80;
t242 = t240+t241-lambda.*t32.*t77-mu.*t28.*t75;
t243 = lambda.*t6.*t77;
t244 = mu.*t12.*t75;
t245 = t243+t244;
t246 = t48.*t242.*(-1.0./6.0)-t51.*t52.*t245.*(1.0./6.0);
t247 = lambda.*t2.*t79;
t248 = mu.*t28.*t82;
t249 = t247+t248;
t250 = t48.*t249.*(1.0./6.0);
t251 = lambda.*t2.*t77;
t252 = t51.*t52.*(t251-mu.*t12.*t82).*(1.0./6.0);
t253 = t250+t252;
t254 = t8.*t64.*t79;
t357 = mu.*t68.*t75;
t255 = t226+t254-t357-t8.*t67.*t77;
t256 = t8.*t64.*t77;
t257 = t192+t228+t256;
t258 = t48.*t255.*(-1.0./6.0)-t51.*t52.*t257.*(1.0./6.0);
t259 = lambda.*t61.*t79;
t260 = mu.*t64.*t80;
t261 = t259+t260-lambda.*t68.*t77-mu.*t67.*t75;
t262 = t48.*t261.*(1.0./6.0);
t263 = lambda.*t61.*t77;
t264 = mu.*t64.*t75;
t265 = t51.*t52.*(t263+t264).*(1.0./6.0);
t266 = t262+t265;
t267 = lambda.*t72.*t79;
t268 = t48.*(t267-mu.*t67.*t82).*(1.0./6.0);
t269 = lambda.*(t10-t34-t63+t70+t76-t81).*(t50-t55+t56-t58+t71-t73);
t270 = mu.*t64.*t82;
t271 = t269+t270;
t272 = t51.*t52.*t271.*(1.0./6.0);
t273 = t268+t272;
t274 = t4-t30-t60+t66+t74-t78;
t275 = t10-t34-t63+t70+t76-t81;
t276 = mu.*t6.*t79;
t277 = lambda.*t12.*t80;
t278 = t276+t277-lambda.*t28.*t75-mu.*t32.*t77;
t279 = lambda.*t12.*t75;
t280 = mu.*t6.*t77;
t281 = t279+t280;
t282 = t48.*t278.*(-1.0./6.0)-t51.*t52.*t281.*(1.0./6.0);
t283 = t6.*t8.*t80;
t341 = mu.*t28.*t77;
t284 = t127+t283-t341-t8.*t32.*t75;
t285 = t48.*t284.*(1.0./6.0);
t287 = t6.*t8.*t75;
t288 = lambda.*t2.*t80;
t289 = mu.*t32.*t82;
t290 = t288+t289;
t291 = lambda.*t2.*t75;
t292 = t291-mu.*t6.*t82;
t293 = t48.*t290.*(-1.0./6.0)-t51.*t52.*t292.*(1.0./6.0);
t294 = lambda.*t67.*t75;
t295 = mu.*t68.*t77;
t296 = t294+t295-lambda.*t64.*t80-mu.*t61.*t79;
t297 = lambda.*t64.*t75;
t298 = mu.*t61.*t77;
t299 = t51.*t52.*(t297+t298).*(1.0./6.0);
t300 = t299-t48.*t296.*(1.0./6.0);
t301 = t8.*t61.*t80;
t358 = mu.*t67.*t77;
t302 = t227+t301-t358-t8.*t68.*t75;
t303 = t8.*t61.*t75;
t304 = t192+t229+t303;
t305 = t48.*t302.*(-1.0./6.0)-t51.*t52.*t304.*(1.0./6.0);
t306 = lambda.*t72.*t80;
t307 = t306-mu.*t68.*t82;
t308 = mu.*t61.*t82;
t309 = lambda.*t72.*t75;
t310 = t308+t309;
t311 = t48.*t307.*(-1.0./6.0)-t51.*t52.*t310.*(1.0./6.0);
t312 = lambda.*t75.*t79;
t313 = mu.*t75.*t79;
t314 = lambda.*t77.*t80;
t315 = mu.*t77.*t80;
t316 = t48.*(t312+t313+t314+t315).*(1.0./6.0);
t317 = lambda.*(t4-t30-t60+t66+t74-t78).*(t10-t34-t63+t70+t76-t81);
t318 = mu.*(t4-t30-t60+t66+t74-t78).*(t10-t34-t63+t70+t76-t81);
t319 = t317+t318;
t320 = t51.*t52.*t319.*(1.0./6.0);
t321 = t316+t320;
t322 = t82.^2;
t323 = mu.*t322;
t324 = t10-t34-t63+t70+t76-t81;
t325 = t4-t30-t60+t66+t74-t78;
t326 = lambda.*t28.*t82;
t327 = mu.*t2.*t79;
t328 = t326+t327;
t329 = t48.*t328.*(1.0./6.0);
t330 = lambda.*t12.*t82;
t331 = t330-mu.*t2.*t77;
t332 = t329-t51.*t52.*t331.*(1.0./6.0);
t333 = lambda.*t32.*t82;
t334 = mu.*t2.*t80;
t335 = t333+t334;
t336 = lambda.*t6.*t82;
t337 = t336-mu.*t2.*t75;
t338 = t51.*t52.*t337.*(1.0./6.0);
t339 = t338-t48.*t335.*(1.0./6.0);
t342 = t51.*t52.*(t128+t129-t2.*t8.*t82).*(1.0./6.0);
t343 = lambda.*t67.*t82;
t344 = t343-mu.*t72.*t79;
t345 = lambda.*t64.*t82;
t346 = mu.*(t10-t34-t63+t70+t76-t81).*(t50-t55+t56-t58+t71-t73);
t347 = t345+t346;
t348 = t51.*t52.*t347.*(1.0./6.0);
t349 = t348-t48.*t344.*(1.0./6.0);
t350 = lambda.*t68.*t82;
t351 = t350-mu.*t72.*t80;
t352 = t48.*t351.*(1.0./6.0);
t353 = lambda.*t61.*t82;
t354 = mu.*t72.*t75;
t355 = t353+t354;
t356 = t352-t51.*t52.*t355.*(1.0./6.0);
t359 = t8.*t72.*t82;
t360 = t228+t229+t359;
t361 = lambda.*t79.*t82;
t362 = mu.*t79.*t82;
t363 = t361+t362;
t364 = t48.*t363.*(1.0./6.0);
t365 = lambda.*t77.*t82;
t366 = mu.*t77.*t82;
t367 = t51.*t52.*(t365+t366).*(1.0./6.0);
t368 = t364+t367;
t369 = lambda.*t80.*t82;
t370 = mu.*t80.*t82;
t371 = t369+t370;
t372 = lambda.*t75.*t82;
t373 = mu.*t75.*t82;
t374 = t372+t373;
t375 = t48.*t371.*(-1.0./6.0)-t51.*t52.*t374.*(1.0./6.0);
t376 = mu.*t75.*t80.*2.0;
t377 = mu.*t77.*t79.*2.0;
t378 = t4-t30-t60+t66+t74-t78;
t379 = t10-t34-t63+t70+t76-t81;
t380 = mu.*t80.*t83;
t381 = mu.*t79.*t84;
t382 = mu.*t75.*t83;
t383 = mu.*t77.*t84;
t384 = t8.*t28.*t84;
t385 = t130+t384;
t386 = t8.*t12.*t84;
t387 = t99+t132+t386;
t388 = t51.*t52.*t387.*(1.0./6.0);
t389 = t388-t48.*t385.*(1.0./6.0);
t390 = mu.*t28.*t83;
t391 = lambda.*t32.*t84;
t392 = t390+t391;
t393 = t48.*t392.*(1.0./6.0);
t394 = lambda.*t6.*t84;
t395 = mu.*t12.*t83;
t396 = t394+t395;
t397 = t393-t51.*t52.*t396.*(1.0./6.0);
t398 = lambda.*t2.*t84;
t399 = mu.*t12.*t51;
t400 = t398+t399;
t401 = t51.*t52.*t400.*(1.0./6.0);
t402 = t401-mu.*t28.*t48.*t51.*(1.0./6.0);
t403 = t8.*t67.*t84;
t404 = t230+t403;
t405 = t48.*t404.*(1.0./6.0);
t406 = t8.*t64.*t84;
t468 = mu.*t51.*t72;
t407 = t232+t406-t468;
t408 = t405-t51.*t52.*t407.*(1.0./6.0);
t409 = mu.*t67.*t83;
t410 = lambda.*t68.*t84;
t411 = t409+t410;
t412 = lambda.*t61.*t84;
t413 = mu.*t64.*t83;
t414 = t412+t413;
t415 = t51.*t52.*t414.*(1.0./6.0);
t416 = t415-t48.*t411.*(1.0./6.0);
t417 = lambda.*t72.*t84;
t418 = t51.*t52.*(t417-mu.*t51.*t64).*(1.0./6.0);
t419 = mu.*t48.*t51.*t67.*(1.0./6.0);
t420 = t418+t419;
t421 = t8.*t79.*t84;
t422 = t380+t421;
t423 = t8.*t77.*t84;
t483 = mu.*t51.*t82;
t424 = t382+t423-t483;
t425 = t48.*t422.*(-1.0./6.0)-t51.*t52.*t424.*(1.0./6.0);
t426 = mu.*t79.*t83;
t427 = lambda.*t80.*t84;
t428 = t426+t427;
t429 = t48.*t428.*(1.0./6.0);
t430 = lambda.*t75.*t84;
t431 = mu.*t77.*t83;
t432 = t51.*t52.*(t430+t431).*(1.0./6.0);
t433 = t429+t432;
t434 = lambda.*t82.*t84;
t435 = t434-mu.*t51.*t77;
t436 = t51.*t52.*t435.*(1.0./6.0);
t437 = t436-mu.*t48.*t51.*t79.*(1.0./6.0);
t438 = lambda.*t28.*t83;
t439 = mu.*t32.*t84;
t440 = t438+t439;
t441 = t48.*t440.*(1.0./6.0);
t442 = lambda.*t12.*t83;
t443 = mu.*t6.*t84;
t444 = t442+t443;
t445 = t441-t51.*t52.*t444.*(1.0./6.0);
t446 = t8.*t32.*t83;
t447 = t131+t446;
t448 = t6.*t8.*t83;
t449 = t99+t133+t448;
t450 = t51.*t52.*t449.*(1.0./6.0);
t451 = t450-t48.*t447.*(1.0./6.0);
t452 = lambda.*t2.*t83;
t453 = mu.*t6.*t51;
t454 = t452+t453;
t455 = mu.*t32.*t48.*t51.*(1.0./6.0);
t456 = t455-t51.*t52.*t454.*(1.0./6.0);
t457 = lambda.*t67.*t83;
t458 = mu.*t68.*t84;
t459 = t457+t458;
t460 = lambda.*t64.*t83;
t461 = mu.*t61.*t84;
t462 = t460+t461;
t463 = t51.*t52.*t462.*(1.0./6.0);
t464 = t463-t48.*t459.*(1.0./6.0);
t465 = t8.*t68.*t83;
t466 = t231+t465;
t467 = t48.*t466.*(1.0./6.0);
t469 = t8.*t61.*t83;
t470 = lambda.*t72.*t83;
t471 = t470-mu.*t51.*t61;
t472 = t51.*t52.*t471.*(-1.0./6.0)-mu.*t48.*t51.*t68.*(1.0./6.0);
t473 = lambda.*t79.*t83;
t474 = mu.*t80.*t84;
t475 = t473+t474;
t476 = t48.*t475.*(1.0./6.0);
t477 = lambda.*t77.*t83;
t478 = mu.*t75.*t84;
t479 = t51.*t52.*(t477+t478).*(1.0./6.0);
t480 = t476+t479;
t481 = t8.*t80.*t83;
t482 = t381+t481;
t484 = t8.*t75.*t83;
t485 = lambda.*t82.*t83;
t486 = t485-mu.*t51.*t75;
t487 = mu.*t48.*t51.*t80.*(1.0./6.0);
t488 = t487-t51.*t52.*t486.*(1.0./6.0);
t489 = lambda.*t83.*t84;
t490 = mu.*t83.*t84;
t491 = t489+t490;
t492 = t51.*t52.*t491.*(1.0./6.0);
t493 = t51.^2;
t494 = mu.*t493;
t495 = t84.^2;
t496 = t83.^2;
t497 = lambda.*t12.*t51;
t498 = mu.*t2.*t84;
t499 = t497+t498;
t500 = t51.*t52.*t499.*(1.0./6.0);
t501 = t500-lambda.*t28.*t48.*t51.*(1.0./6.0);
t502 = lambda.*t6.*t51;
t503 = mu.*t2.*t83;
t504 = t502+t503;
t505 = lambda.*t32.*t48.*t51.*(1.0./6.0);
t506 = t505-t51.*t52.*t504.*(1.0./6.0);
t507 = t130+t131;
t508 = t2.*t8.*t51;
t509 = t132+t133+t508;
t510 = t51.*t52.*t509.*(1.0./6.0);
t511 = t510-t48.*t507.*(1.0./6.0);
t512 = lambda.*t51.*t64;
t513 = t512-mu.*t72.*t84;
t514 = lambda.*t48.*t51.*t67.*(1.0./6.0);
t515 = t514-t51.*t52.*t513.*(1.0./6.0);
t516 = lambda.*t51.*t61;
t517 = t516-mu.*t72.*t83;
t518 = t51.*t52.*t517.*(1.0./6.0);
t519 = t518-lambda.*t48.*t51.*t68.*(1.0./6.0);
t520 = t230+t231;
t521 = t48.*t520.*(1.0./6.0);
t522 = t232+t233-t8.*t51.*t72;
t523 = t521-t51.*t52.*t522.*(1.0./6.0);
t524 = lambda.*t51.*t77;
t525 = t524-mu.*t82.*t84;
t526 = t51.*t52.*t525.*(-1.0./6.0)-lambda.*t48.*t51.*t79.*(1.0./6.0);
t527 = lambda.*t51.*t75;
t528 = t51.*t52.*(t527-mu.*t82.*t83).*(1.0./6.0);
t529 = lambda.*t48.*t51.*t80.*(1.0./6.0);
t530 = t528+t529;
t531 = t380+t381;
t532 = t382+t383-t8.*t51.*t82;
t533 = t48.*t531.*(-1.0./6.0)-t51.*t52.*t532.*(1.0./6.0);
t534 = lambda.*t51.*t84;
t535 = mu.*t51.*t84;
t536 = t534+t535;
t537 = lambda.*t51.*t83;
t538 = mu.*t51.*t83;
t539 = t537+t538;
t540 = t51.*t52.*t539.*(1.0./6.0);
t541 = mu.*t496;
t542 = mu.*t495;
der_st_dz4 = reshape([t48.*(t116+t8.*t12.*t28.*2.0).*(1.0./6.0)-t51.*t52.*(t96+t118+t8.*t97).*(1.0./6.0),t94,t107,t140,t167,t197,t239,t282,t332,t389,t445,t501,t94,t48.*(t117+t6.*t8.*t32.*2.0).*(1.0./6.0)-t51.*t52.*(t96+t119+t8.*t98).*(1.0./6.0),t115,t150,t48.*t170.*(-1.0./6.0)+t51.*t52.*(t125-t171+t172).*(1.0./6.0),t202,t246,t285+t51.*t52.*(t129-t286+t287).*(1.0./6.0),t339,t397,t451,t506,t107,t115,t48.*(t116+t117).*(1.0./6.0)-t51.*t52.*(t118+t119+t8.*t95).*(1.0./6.0),t156,t178,t206,t253,t293,t342+t48.*(t126+t127-t340-t341).*(1.0./6.0),t402,t456,t511,t140,t150,t156,t48.*(t221+t8.*t64.*t67.*2.0).*(1.0./6.0)-t51.*t52.*(t223+t8.*t190+mu.*t157.^2).*(1.0./6.0),t188,t213,t258,t300,t349,t408,t464,t515,t167,t48.*t170.*(-1.0./6.0)+t51.*t52.*(t125+t172-mu.*t2.*t72).*(1.0./6.0),t178,t188,t48.*(t222+t8.*t61.*t68.*2.0).*(1.0./6.0)-t51.*t52.*(t224+t8.*t191+mu.*t189.^2).*(1.0./6.0),t220,t266,t305,t356,t416,t467-t51.*t52.*(t233-t468+t469).*(1.0./6.0),t519,t197,t202,t206,t213,t220,t48.*(t221+t222).*(1.0./6.0)-t51.*t52.*(t223+t224+t8.*t225.^2).*(1.0./6.0),t273,t311,t48.*(t226+t227-t357-t358).*(-1.0./6.0)-t51.*t52.*t360.*(1.0./6.0),t420,t472,t523,t239,t246,t253,t258,t266,t273,t48.*(t376+t8.*t77.*t79.*2.0).*(-1.0./6.0)-t51.*t52.*(t323+mu.*t274.^2+t8.*t275.^2).*(1.0./6.0),t321,t368,t425,t480,t526,t282,t285+t51.*t52.*(t129+t287-mu.*t2.*t82).*(1.0./6.0),t293,t300,t305,t311,t321,t48.*(t377+t8.*t75.*t80.*2.0).*(-1.0./6.0)-t51.*t52.*(t323+mu.*t324.^2+t8.*t325.^2).*(1.0./6.0),t375,t433,t48.*t482.*(-1.0./6.0)-t51.*t52.*(t383-t483+t484).*(1.0./6.0),t530,t332,t339,t342+t48.*(t126+t127-mu.*t28.*t77-mu.*t32.*t75).*(1.0./6.0),t349,t356,t48.*(t226+t227-mu.*t68.*t75-mu.*t67.*t77).*(-1.0./6.0)-t51.*t52.*t360.*(1.0./6.0),t368,t375,t48.*(t376+t377).*(-1.0./6.0)-t51.*t52.*(t8.*t322+mu.*t378.^2+mu.*t379.^2).*(1.0./6.0),t437,t488,t533,t389,t397,t402,t408,t416,t420,t425,t433,t437,t51.*t52.*(t494+t541+t8.*t495).*(-1.0./6.0),t492,t51.*t52.*t536.*(-1.0./6.0),t445,t451,t456,t464,t467-t51.*t52.*(t233+t469-mu.*t51.*t72).*(1.0./6.0),t472,t480,t48.*t482.*(-1.0./6.0)-t51.*t52.*(t383+t484-mu.*t51.*t82).*(1.0./6.0),t488,t492,t51.*t52.*(t494+t542+t8.*t496).*(-1.0./6.0),t540,t501,t506,t511,t515,t519,t523,t526,t530,t533,t51.*t52.*t536.*(-1.0./6.0),t540,t51.*t52.*(t541+t542+t8.*t493).*(-1.0./6.0)],[12, 12]);
