// Benchmark "ex" written by ABC on Fri Jul  4 15:14:12 2025

module top (
    a, b, c, d,
    F0  );
  input  a, b, c, d;
  output F0;
  assign F0 = ~c & ~d;
endmodule


