clear all
close all
clc

data={
  % abc sequences
  'abcabcabcabcabcabcabcabcabcabcabcabc'
  'bcabcabcabcabcabcabcabcabcabcabcabc'
  'bcabcabcabcabcabcabcabcabcabcabcabc'
  'cabcabcabcabcabcabcabcabcabcabcabc'
  'abcabcabc'
  'bcabcabc'
  'abcabcabcabc'
  'bcabcabcabcabcab' 
  'abcabcabcabcabcabc'
  % acb sequences
  'acbacbacbacbacbacbacbacb'
  'bacbacbacbacbacbac'
  'acbacbacbabcbacbacbacbacbacbacbacbacbac'
  'acbacbacbabcbacbacbacbacbacbacbacbacbac'
  'cbacbabcbacbacbacbacbacbacbacbacbac'
  'acbabcbacbacbacbacbacbacbacbacbac'
  'acbacbacb'
  % random sequences
  'aabbbabbbababbbaabbbbababbbaaa'
  'aabbababbababbbbabbaaababaabaa'
  'abbabbbbbbbaabbabbaaaaaabababa'
  'baabaabbabaaaabbabaaabbaabbbaa'
  'bbaaabbababaababbbbbaaabaaabba'
  };

disp('Showing data...');
data
disp(['The data is made of abc sequences, acb sequences, and random ab' ...
      ' sequences']);
data1 = data(1);
net = vbhmm(data1,'abc',2,100,1e-6); % 2 states
Fv = vbhmm_cF(data,'abc',net);

