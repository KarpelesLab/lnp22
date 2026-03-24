package modmath

import "testing"

const testQ int64 = 8380417

func TestReduce(t *testing.T) {
	tests := []struct {
		a, q, want int64
	}{
		{0, testQ, 0},
		{1, testQ, 1},
		{-1, testQ, testQ - 1},
		{testQ, testQ, 0},
		{testQ + 1, testQ, 1},
		{-testQ, testQ, 0},
		{2 * testQ, testQ, 0},
		{-2*testQ - 3, testQ, testQ - 3},
	}
	for _, tt := range tests {
		got := Reduce(tt.a, tt.q)
		if got != tt.want {
			t.Errorf("Reduce(%d, %d) = %d, want %d", tt.a, tt.q, got, tt.want)
		}
	}
}

func TestCenterReduce(t *testing.T) {
	half := testQ / 2
	tests := []struct {
		a, q, want int64
	}{
		{0, testQ, 0},
		{1, testQ, 1},
		{half, testQ, half},
		{half + 1, testQ, -(half)},
		{-1, testQ, -1},
		{testQ, testQ, 0},
	}
	for _, tt := range tests {
		got := CenterReduce(tt.a, tt.q)
		if got != tt.want {
			t.Errorf("CenterReduce(%d, %d) = %d, want %d", tt.a, tt.q, got, tt.want)
		}
	}
}

func TestModMul(t *testing.T) {
	tests := []struct {
		a, b, q, want int64
	}{
		{0, 5, testQ, 0},
		{1, 5, testQ, 5},
		{2, 3, testQ, 6},
		{testQ - 1, testQ - 1, testQ, 1}, // (-1)*(-1) = 1
		{1000000, 1000000, testQ, Reduce(1000000*1000000, testQ)},
	}
	for _, tt := range tests {
		got := ModMul(tt.a, tt.b, tt.q)
		if got != tt.want {
			t.Errorf("ModMul(%d, %d, %d) = %d, want %d", tt.a, tt.b, tt.q, got, tt.want)
		}
	}
}

func TestModExp(t *testing.T) {
	tests := []struct {
		base, exp, q, want int64
	}{
		{2, 0, testQ, 1},
		{2, 1, testQ, 2},
		{2, 10, testQ, 1024},
		{3, testQ - 1, testQ, 1}, // Fermat's little theorem
		{7, testQ - 1, testQ, 1},
	}
	for _, tt := range tests {
		got := ModExp(tt.base, tt.exp, tt.q)
		if got != tt.want {
			t.Errorf("ModExp(%d, %d, %d) = %d, want %d", tt.base, tt.exp, tt.q, got, tt.want)
		}
	}
}

func TestModInverse(t *testing.T) {
	vals := []int64{1, 2, 3, 7, 1000, testQ - 1, 123456}
	for _, a := range vals {
		inv := ModInverse(a, testQ)
		prod := ModMul(a, inv, testQ)
		if prod != 1 {
			t.Errorf("ModInverse(%d, %d) = %d, but %d*%d mod %d = %d",
				a, testQ, inv, a, inv, testQ, prod)
		}
	}
}

func TestModInversePanicsOnZero(t *testing.T) {
	defer func() {
		if r := recover(); r == nil {
			t.Error("ModInverse(0, q) should panic")
		}
	}()
	ModInverse(0, testQ)
}

func TestAbs(t *testing.T) {
	tests := []struct {
		a, want int64
	}{
		{0, 0}, {5, 5}, {-5, 5}, {-1, 1},
	}
	for _, tt := range tests {
		if got := Abs(tt.a); got != tt.want {
			t.Errorf("Abs(%d) = %d, want %d", tt.a, got, tt.want)
		}
	}
}
