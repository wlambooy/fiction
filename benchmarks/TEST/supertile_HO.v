module top (
    A,     // Input for FO1 and FO2
    B,     // Input for FO1 and FO2
    AND_out,  // Output from AND gate
    OR_out    // Output from OR gate
);

input A, B;
output AND_out, OR_out;

// FO1 Logic
wire FO1_buf;  // Buffered signal from FO1
wire FO1_ho;   // Direct signal from FO1

// FO2 Logic
wire FO2_ho;  // Direct signal from FO2
wire FO2_inv;   // Inverted signal from FO2

// FO1 logic
assign FO1_buf = A;   // Buffered input A for FO1 (BUF)
assign FO1_ho = A;     // Direct connection from FO1

// FO2 logic
assign FO2_ho = B;    // Direct connection from FO2
assign FO2_inv = ~B;  // Inverted input B for FO2 (INV)

// AND gate logic
assign AND_out = FO1_buf & FO1_ho;  // FO1 → BUF → AND, FO1 → HO → AND

// OR gate logic
assign OR_out = FO2_ho | FO2_inv;   // FO2 → HO → OR, FO2 → INV → OR

endmodule