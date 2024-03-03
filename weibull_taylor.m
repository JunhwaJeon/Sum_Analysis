close all;
clear all;
clc;

syms kap mu omeg t
fun=exp(kap*mu*t)*(t^mu)*(1-t)^(-2);
fun_int=vpaintegral(fun,t,0,1);
fun_int