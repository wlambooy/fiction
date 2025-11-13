// Benchmark "ex" written by ABC on Fri Jul  4 15:14:11 2025

module top (
    a, b, c,
    F0  );
  input  a, b, c;
  output F0;
  wire new_n5, new_n6, new_n7, new_n8, new_n9, new_n10, new_n11, new_n12;
  assign new_n5 = ~b & c;
  assign new_n6 = b & ~c;
  assign new_n7 = ~new_n5 & ~new_n6;
  assign new_n8 = a & ~new_n7;
  assign new_n9 = b & c;
  assign new_n10 = ~b & ~c;
  assign new_n11 = ~new_n9 & ~new_n10;
  assign new_n12 = ~a & ~new_n11;
  assign F0 = new_n8 | new_n12;
endmodule


