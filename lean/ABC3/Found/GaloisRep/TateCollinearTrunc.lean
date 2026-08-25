import ABC3.Found.GaloisRep.TateCollinearAnalytic

/-!
# Galois (G6) 第 250 ブロック —— **★★★★★★★★共線性を有限の主張の族に落とす**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★第 222–223 の骨格をそのまま 3 点に使う

葉 (b)(Tate 方程式)では

    差 → 切り詰め → 万有な環での整除性

という 3 段で有限化した(第 222・223)。共線性でも**同じ骨格が使える**。
違うのは点が 3 つになることだけである。

3 点は `(u, vw)`, `(v, uw)`, `(w, uv)` で、`q = u v w`。★各点の相方が多項式になるので
**追加の変数は要らない**(§9-562)。

    collDefect u v w q := X₁(Y₂ − Y₃) + X₂(Y₃ − Y₁) + X₃(Y₁ − Y₂)

## ★★★★行列式なので互換で符号が変わる

    collDefectTrunc n v u w q = −collDefectTrunc n u v w q     (`collDefectTrunc_swap12`)
    collDefectTrunc n u w v q = −collDefectTrunc n u v w q     (`collDefectTrunc_swap23`)

★これは `ring` で出る(`v * u = u * v` を先に潰すだけ)。第 224 の `tateDefectTrunc_symm`
が `A` 側の整除性を無料にしたのと同じで、ここでも **`U` 側を示せば `V`・`W` 側は
対称性で済む**。

## ★★切り詰めの函手性は在庫を使う

`map_tateXtrunc`・`map_tateYtrunc`(第 222)を 3 点それぞれに当てるだけ。
単元条件は 12 個になるので `CollUnits` にまとめた。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `collDefect` | ★★★★★★共線性の差(形式側) |
| `collDefectTrunc` | ★★★その `n` 次の切り詰め |
| `collDefect_eq_zero_of_trunc` | ★★★★★★★有限の主張の族に落ちる |
| `collDefectTrunc_swap12`・`_swap23` | ★★★★互換で符号が変わる |
| `map_collDefectTrunc` | ★★★★★★環準同型で移る |
| `collDefect_eq_zero_of_specialize` | ★★★★★★★★**特殊化の規準** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★★★共線性の差と切り詰め -/

/-- ★★★★★★**共線性の差(形式側)**——3 点は `(u, vw)`, `(v, uw)`, `(w, uv)`、`q = uvw`。 -/
noncomputable def collDefect [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I) : R :=
  tateXpair u (v * w) q hq * (tateYpair v (u * w) q hq - tateYpair w (u * v) q hq)
    + tateXpair v (u * w) q hq * (tateYpair w (u * v) q hq - tateYpair u (v * w) q hq)
    + tateXpair w (u * v) q hq * (tateYpair u (v * w) q hq - tateYpair v (u * w) q hq)

/-- ★★★**共線性の差の `n` 次の切り詰め**——無限和も極限も含まない。 -/
noncomputable def collDefectTrunc (n : ℕ) (u v w q : R) : R :=
  tateXtrunc n u (v * w) q * (tateYtrunc n v (u * w) q - tateYtrunc n w (u * v) q)
    + tateXtrunc n v (u * w) q * (tateYtrunc n w (u * v) q - tateYtrunc n u (v * w) q)
    + tateXtrunc n w (u * v) q * (tateYtrunc n u (v * w) q - tateYtrunc n v (u * w) q)

theorem collDefect_sub_trunc [IsAdicComplete I R] (n : ℕ) (u v w q : R) (hq : q ∈ I) :
    collDefect u v w q hq - collDefectTrunc n u v w q ∈ I ^ n := by
  refine Ideal.Quotient.eq.1 ?_
  have hX1 : (Ideal.Quotient.mk (I ^ n)) (tateXpair u (v * w) q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateXtrunc n u (v * w) q) :=
    Ideal.Quotient.eq.2 (tateXpair_sub_trunc n u (v * w) q hq)
  have hX2 : (Ideal.Quotient.mk (I ^ n)) (tateXpair v (u * w) q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateXtrunc n v (u * w) q) :=
    Ideal.Quotient.eq.2 (tateXpair_sub_trunc n v (u * w) q hq)
  have hX3 : (Ideal.Quotient.mk (I ^ n)) (tateXpair w (u * v) q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateXtrunc n w (u * v) q) :=
    Ideal.Quotient.eq.2 (tateXpair_sub_trunc n w (u * v) q hq)
  have hY1 : (Ideal.Quotient.mk (I ^ n)) (tateYpair u (v * w) q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateYtrunc n u (v * w) q) :=
    Ideal.Quotient.eq.2 (tateYpair_sub_trunc n u (v * w) q hq)
  have hY2 : (Ideal.Quotient.mk (I ^ n)) (tateYpair v (u * w) q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateYtrunc n v (u * w) q) :=
    Ideal.Quotient.eq.2 (tateYpair_sub_trunc n v (u * w) q hq)
  have hY3 : (Ideal.Quotient.mk (I ^ n)) (tateYpair w (u * v) q hq)
      = (Ideal.Quotient.mk (I ^ n)) (tateYtrunc n w (u * v) q) :=
    Ideal.Quotient.eq.2 (tateYpair_sub_trunc n w (u * v) q hq)
  simp only [collDefect, collDefectTrunc, map_sub, map_add, map_mul,
    hX1, hX2, hX3, hY1, hY2, hY3]

/-- ★★★★★★★**共線性も有限の代数的主張の族に落ちる**。 -/
theorem collDefect_eq_zero_of_trunc [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (h : ∀ n, collDefectTrunc n u v w q ∈ I ^ n) : collDefect u v w q hq = 0 := by
  refine eq_zero_of_mem_pow (I := I) fun n => ?_
  have hd := collDefect_sub_trunc (I := I) n u v w q hq
  have hsum := Ideal.add_mem (I ^ n) hd (h n)
  simpa using hsum

/-! ## ★★★★互換で符号が変わる -/

/-- ★★★★**行列式は互換で符号が変わる**(1 と 2 の入れ替え)。

★これで `U` 側の整除性から `V` 側が無料で出る。 -/
theorem collDefectTrunc_swap12 (n : ℕ) (u v w q : R) :
    collDefectTrunc n v u w q = -collDefectTrunc n u v w q := by
  simp only [collDefectTrunc, show v * u = u * v from mul_comm v u]
  ring

/-- ★★★★**行列式は互換で符号が変わる**(2 と 3 の入れ替え)。 -/
theorem collDefectTrunc_swap23 (n : ℕ) (u v w q : R) :
    collDefectTrunc n u w v q = -collDefectTrunc n u v w q := by
  simp only [collDefectTrunc, show w * v = v * w from mul_comm w v]
  ring

end ABC3.Found.GaloisRep

/-! ## ★★★★★★切り詰めの函手性(完備性は要らない) -/

namespace ABC3.Found.GaloisRep

section Functorial

variable {R R' : Type} [CommRing R] [CommRing R']

/-- ★★★★★★**切り詰めた共線性の差は任意の環準同型で移る**——第 222 の
`map_tateXtrunc`・`map_tateYtrunc` を 3 点それぞれに当てるだけ。 -/
theorem map_collDefectTrunc (φ : R →+* R') (n : ℕ) (u v w q : R)
    (hu : IsUnit (1 - u)) (hv : IsUnit (1 - v)) (hw : IsUnit (1 - w))
    (hvw : IsUnit (1 - v * w)) (huw : IsUnit (1 - u * w)) (huv : IsUnit (1 - u * v))
    (hqu : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * u))
    (hqv : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * v))
    (hqw : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * w))
    (hqvw : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * (v * w)))
    (hquw : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * (u * w)))
    (hquv : ∀ m, m < n → IsUnit (1 - q ^ (m + 1) * (u * v))) :
    φ (collDefectTrunc n u v w q) = collDefectTrunc n (φ u) (φ v) (φ w) (φ q) := by
  simp only [collDefectTrunc, map_sub, map_add, map_mul,
    map_tateXtrunc φ n u (v * w) q hu hvw hqu hqvw,
    map_tateXtrunc φ n v (u * w) q hv huw hqv hquw,
    map_tateXtrunc φ n w (u * v) q hw huv hqw hquv,
    map_tateYtrunc φ n u (v * w) q hu hvw hqu hqvw,
    map_tateYtrunc φ n v (u * w) q hv huw hqv hquw,
    map_tateYtrunc φ n w (u * v) q hw huv hqw hquv]

end Functorial

/-! ## ★★★★★★★★特殊化の規準 -/

/-- ★3 点すべての切り詰めが多項式として意味を持つ条件(12 個)。 -/
structure CollUnits {S : Type} [CommRing S] (U V W : S) : Prop where
  hu : IsUnit (1 - U)
  hv : IsUnit (1 - V)
  hw : IsUnit (1 - W)
  hvw : IsUnit (1 - V * W)
  huw : IsUnit (1 - U * W)
  huv : IsUnit (1 - U * V)
  hqu : ∀ m : ℕ, IsUnit (1 - (U * V * W) ^ (m + 1) * U)
  hqv : ∀ m : ℕ, IsUnit (1 - (U * V * W) ^ (m + 1) * V)
  hqw : ∀ m : ℕ, IsUnit (1 - (U * V * W) ^ (m + 1) * W)
  hqvw : ∀ m : ℕ, IsUnit (1 - (U * V * W) ^ (m + 1) * (V * W))
  hquw : ∀ m : ℕ, IsUnit (1 - (U * V * W) ^ (m + 1) * (U * W))
  hquv : ∀ m : ℕ, IsUnit (1 - (U * V * W) ^ (m + 1) * (U * V))

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★★★**特殊化の規準(共線性)**——万有な環 `S` の上で切り詰めが
`(UVW)^n` の倍数なら、任意の完備環で 3 点は共線になる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem collDefect_eq_zero_of_specialize [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q)
    {S : Type} [CommRing S] (U V W : S) (φ : S →+* R)
    (hU : φ U = u) (hV : φ V = v) (hW : φ W = w)
    (hcu : CollUnits U V W)
    (H : ∀ n : ℕ, ((U * V * W) ^ n) ∣ collDefectTrunc n U V W (U * V * W)) :
    collDefect u v w q hq = 0 := by
  refine collDefect_eq_zero_of_trunc (I := I) u v w q hq fun n => ?_
  have hQ : φ (U * V * W) = q := by rw [map_mul, map_mul, hU, hV, hW, huvw]
  have hmap := map_collDefectTrunc φ n U V W (U * V * W) hcu.hu hcu.hv hcu.hw
    hcu.hvw hcu.huw hcu.huv
    (fun m _ => hcu.hqu m) (fun m _ => hcu.hqv m) (fun m _ => hcu.hqw m)
    (fun m _ => hcu.hqvw m) (fun m _ => hcu.hquw m) (fun m _ => hcu.hquv m)
  rw [hU, hV, hW, hQ] at hmap
  obtain ⟨c, hc⟩ := H n
  rw [← hmap, hc, map_mul, map_pow, hQ]
  exact Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)

/-! ## ★出典の紐付け(`.src`) -/

def collDefect.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——共線性の差)",
    sectionId := "genell-def-3-3" }

def collDefect_eq_zero_of_specialize.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——共線性の特殊化の規準)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
