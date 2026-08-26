import ABC3.Found.GaloisRep.TateIdentify

/-!
# Galois (G6) 第 265 ブロック —— **★★★★★★★葉 (e) の縮小性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★誤差項の差は 1 つ位が上がる

葉 (e) の存在の段は `Λ(a) = c` の解 `a` を作ることである。第 263 より
`Λ(a) = a + κ(a)`(`κ(a) ∈ I`)なので、`a = c − κ(a)` の不動点を探せばよい。
★そのために要るのは **`κ` の変動が 1 つ位を上げる**こと、すなわち

    a − b ∈ I^j  ⟹  ε(a) − ε(b) ∈ I^{j+1},   ε(a) := Y − a(X+Y)

である(`tateEps_diff_mem`)。

## ★★★相方 `w = q/a` は 1 つ位を稼ぐ

`a` が単元なので `w := q·a⁻¹` は `R` の中で取れる(`wOf`)。逆元の差は

    inv(a) − inv(b) = inv(a)·inv(b)·(b − a)  ∈ I^j

なので、`q ∈ I` を掛けて **`w(a) − w(b) ∈ I^{j+1}`**(`wOf_diff_mem`)。
★★これが縮小性の源である——`w` 側の項は最初から 1 つ位を稼いでいる。

## ★★★★★誤差項の展開

第 262 の `tateYterm_eq_mul`(`g(a) = a(f(a)+g(a))`)で主要項が消えるので

    ε = [Tg(a) − f(w) − Tf(w) − g(w) − Tg(w) + s₁]
        − a·[Tf(a) + Tg(a) − g(w) − Tg(w) − s₁]

★**尾と `w` の項だけ**が残る。`a` そのものへの依存は尾 `Tf(a), Tg(a)` を通してのみで、
それらの差は `q ∈ I` を 1 つ含むので `I^{j+1}`(第 259 の `tateXtail_diff_mem`)。
`a·[…]` の項は `(a−b)·[…]`(`[…] ∈ I`)と `b·([…]−[…]')` に分ける。

## ★単元版の差の補題

第 259 の差の補題は `a, b ∈ I` を仮定していた(`1 − a` が単元であるため)。
`a` が**単元**の場合も要るので、仮定を `IsUnit (1 − a)` に直した版を置いた
(`tateXterm_diff_mem_gen`・`tateYterm_diff_mem_gen`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `inverse_diff_mem` | ★★★★逆元の差 |
| `tateXterm_diff_mem_gen`・`tateYterm_diff_mem_gen` | ★★★★単元版の差 |
| `wOf`・`wOf_diff_mem` | ★★★相方 `q/a` と 1 つ位を稼ぐこと |
| `tateEps_eq` | ★★★★★誤差項の展開 |
| `tateEps_diff_mem` | ★★★★★★★**誤差項の差は 1 つ位が上がる** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★逆元と単元版の差 -/

/-- ★★★★**逆元の差**——`inv(a) − inv(b) = inv(a)inv(b)(b − a)`。 -/
theorem inverse_diff_mem {j : ℕ} {a b : R} (ha : IsUnit a) (hb : IsUnit b)
    (hab : a - b ∈ I ^ j) : Ring.inverse a - Ring.inverse b ∈ I ^ j := by
  have hia : Ring.inverse a * a = 1 := Ring.inverse_mul_cancel _ ha
  have hib : Ring.inverse b * b = 1 := Ring.inverse_mul_cancel _ hb
  have hkey : Ring.inverse a - Ring.inverse b
      = (Ring.inverse a * Ring.inverse b) * (b - a) := by
    calc Ring.inverse a - Ring.inverse b
        = Ring.inverse a * (Ring.inverse b * b) - (Ring.inverse a * a) * Ring.inverse b := by
          rw [hia, hib]; ring
      _ = (Ring.inverse a * Ring.inverse b) * (b - a) := by ring
  rw [hkey]
  refine Ideal.mul_mem_left _ _ ?_
  have h := neg_mem hab
  rwa [neg_sub] at h

/-- ★★★★**`f` の差(単元版)**——仮定を `IsUnit (1 − a)` にした形。 -/
theorem tateXterm_diff_mem_gen {j : ℕ} {a b : R} (ha : IsUnit (1 - a)) (hb : IsUnit (1 - b))
    (hab : a - b ∈ I ^ j) : tateXterm a - tateXterm b ∈ I ^ j := by
  set D : R := (1 - a) ^ 2 * (1 - b) ^ 2 with hDdef
  have huD : IsUnit D := by
    simp only [hDdef]; exact (ha.pow 2).mul (hb.pow 2)
  have h1 : Ring.inverse D * D = 1 := Ring.inverse_mul_cancel _ huD
  have hD : D * (tateXterm a - tateXterm b) = (a - b) * (1 - a * b) := by
    simp only [hDdef]
    linear_combination (1 - b) ^ 2 * mul_tateXterm' ha - (1 - a) ^ 2 * mul_tateXterm' hb
  have hkey : tateXterm a - tateXterm b = (a - b) * (Ring.inverse D * (1 - a * b)) := by
    calc tateXterm a - tateXterm b
        = Ring.inverse D * (D * (tateXterm a - tateXterm b)) := by
          rw [← mul_assoc, h1, one_mul]
      _ = Ring.inverse D * ((a - b) * (1 - a * b)) := by rw [hD]
      _ = (a - b) * (Ring.inverse D * (1 - a * b)) := by ring
  rw [hkey]
  exact Ideal.mul_mem_right _ _ hab

/-- ★★★★**`g` の差(単元版)**。 -/
theorem tateYterm_diff_mem_gen {j : ℕ} {a b : R} (ha : IsUnit (1 - a)) (hb : IsUnit (1 - b))
    (hab : a - b ∈ I ^ j) : tateYterm a - tateYterm b ∈ I ^ j := by
  set D : R := (1 - a) ^ 3 * (1 - b) ^ 3 with hDdef
  have huD : IsUnit D := by
    simp only [hDdef]; exact (ha.pow 3).mul (hb.pow 3)
  have h1 : Ring.inverse D * D = 1 := Ring.inverse_mul_cancel _ huD
  have hD : D * (tateYterm a - tateYterm b)
      = (a - b) * (a + b - 3 * (a * b) + a ^ 2 * b ^ 2) := by
    simp only [hDdef]
    linear_combination (1 - b) ^ 3 * mul_tateYterm' ha - (1 - a) ^ 3 * mul_tateYterm' hb
  have hkey : tateYterm a - tateYterm b
      = (a - b) * (Ring.inverse D * (a + b - 3 * (a * b) + a ^ 2 * b ^ 2)) := by
    calc tateYterm a - tateYterm b
        = Ring.inverse D * (D * (tateYterm a - tateYterm b)) := by
          rw [← mul_assoc, h1, one_mul]
      _ = Ring.inverse D * ((a - b) * (a + b - 3 * (a * b) + a ^ 2 * b ^ 2)) := by rw [hD]
      _ = (a - b) * (Ring.inverse D * (a + b - 3 * (a * b) + a ^ 2 * b ^ 2)) := by ring
  rw [hkey]
  exact Ideal.mul_mem_right _ _ hab

/-! ## ★★★相方 `q/a` -/

/-- ★★★`a` が単元のときの相方 `w = q/a`。 -/
noncomputable def wOf (q a : R) : R := q * Ring.inverse a

theorem wOf_mem {q : R} (a : R) (hq : q ∈ I) : wOf q a ∈ I := Ideal.mul_mem_right _ _ hq

theorem mul_wOf (q a : R) (ha : IsUnit a) : a * wOf q a = q := by
  rw [wOf, show a * (q * Ring.inverse a) = q * (a * Ring.inverse a) by ring,
    Ring.mul_inverse_cancel _ ha, mul_one]

/-- ★★★★**相方は 1 つ位を稼ぐ**——`q ∈ I` を掛けるから。 -/
theorem wOf_diff_mem {j : ℕ} {a b q : R} (hq : q ∈ I) (ha : IsUnit a) (hb : IsUnit b)
    (hab : a - b ∈ I ^ j) : wOf q a - wOf q b ∈ I ^ (j + 1) := by
  have h : wOf q a - wOf q b = (Ring.inverse a - Ring.inverse b) * q := by
    rw [wOf, wOf]; ring
  rw [h, pow_succ]
  exact Ideal.mul_mem_mul (inverse_diff_mem ha hb hab) hq

/-! ## ★★★★★★★誤差項 -/

/-- ★★★★★**誤差項 `ε = Y − a(X+Y)` の展開**——主要項が消えて尾と `w` だけ残る。 -/
theorem tateEps_eq [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (hu : IsUnit (1 - a)) :
    tateYpair a w q hq - a * (tateXpair a w q hq + tateYpair a w q hq)
      = (tateYtail a q hq - tateXterm w - tateXtail w q hq - tateYterm w - tateYtail w q hq
          + evalAdic (sigmaSeries 1) q hq)
        - a * (tateXtail a q hq + tateYtail a q hq - tateYterm w - tateYtail w q hq
          - evalAdic (sigmaSeries 1) q hq) := by
  rw [tateXpair, tateYpair]
  linear_combination tateYterm_eq_mul hu

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**誤差項の差は 1 つ位が上がる**——縮小写像の核心。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateEps_diff_mem [IsAdicComplete I R] {j : ℕ} (a b q : R) (hq : q ∈ I)
    (hau : IsUnit a) (hbu : IsUnit b) (ha1 : IsUnit (1 - a)) (hb1 : IsUnit (1 - b))
    (hab : a - b ∈ I ^ j) :
    (tateYpair a (wOf q a) q hq
        - a * (tateXpair a (wOf q a) q hq + tateYpair a (wOf q a) q hq))
      - (tateYpair b (wOf q b) q hq
        - b * (tateXpair b (wOf q b) q hq + tateYpair b (wOf q b) q hq))
      ∈ I ^ (j + 1) := by
  have hwa : wOf q a ∈ I := wOf_mem a hq
  have hwb : wOf q b ∈ I := wOf_mem b hq
  have hwab : wOf q a - wOf q b ∈ I ^ (j + 1) := wOf_diff_mem hq hau hbu hab
  have hq1 : q ∈ I ^ 1 := by rwa [pow_one]
  have hup : ∀ x : R, x ∈ I ^ (j + 1 + 1) → x ∈ I ^ (j + 1) := fun x hx =>
    Ideal.pow_le_pow_right (by omega) hx
  have dTga : tateYtail a q hq - tateYtail b q hq ∈ I ^ (j + 1) :=
    tateYtail_diff_mem (m := 1) a b q hq hq1 hab
  have dTfa : tateXtail a q hq - tateXtail b q hq ∈ I ^ (j + 1) :=
    tateXtail_diff_mem (m := 1) a b q hq hq1 hab
  have dfw : tateXterm (wOf q a) - tateXterm (wOf q b) ∈ I ^ (j + 1) :=
    tateXterm_diff_mem_gen (isUnit_one_sub hwa) (isUnit_one_sub hwb) hwab
  have dgw : tateYterm (wOf q a) - tateYterm (wOf q b) ∈ I ^ (j + 1) :=
    tateYterm_diff_mem_gen (isUnit_one_sub hwa) (isUnit_one_sub hwb) hwab
  have dTfw : tateXtail (wOf q a) q hq - tateXtail (wOf q b) q hq ∈ I ^ (j + 1) :=
    hup _ (tateXtail_diff_mem (m := 1) (wOf q a) (wOf q b) q hq hq1 hwab)
  have dTgw : tateYtail (wOf q a) q hq - tateYtail (wOf q b) q hq ∈ I ^ (j + 1) :=
    hup _ (tateYtail_diff_mem (m := 1) (wOf q a) (wOf q b) q hq hq1 hwab)
  have hE2a : tateXtail a q hq + tateYtail a q hq - tateYterm (wOf q a)
      - tateYtail (wOf q a) q hq - evalAdic (sigmaSeries 1) q hq ∈ I :=
    Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _
      (Ideal.add_mem _ (tateXtail_mem a q hq) (tateYtail_mem a q hq))
      (tateYterm_mem_one hwa)) (tateYtail_mem (wOf q a) q hq))
      (evalAdic_mem (sigmaSeries 1) (by simp [sigmaSeries]) q hq)
  have dE2 : (tateXtail a q hq + tateYtail a q hq - tateYterm (wOf q a)
        - tateYtail (wOf q a) q hq - evalAdic (sigmaSeries 1) q hq)
      - (tateXtail b q hq + tateYtail b q hq - tateYterm (wOf q b)
        - tateYtail (wOf q b) q hq - evalAdic (sigmaSeries 1) q hq) ∈ I ^ (j + 1) := by
    have hkey : (tateXtail a q hq + tateYtail a q hq - tateYterm (wOf q a)
        - tateYtail (wOf q a) q hq - evalAdic (sigmaSeries 1) q hq)
      - (tateXtail b q hq + tateYtail b q hq - tateYterm (wOf q b)
        - tateYtail (wOf q b) q hq - evalAdic (sigmaSeries 1) q hq)
      = (tateXtail a q hq - tateXtail b q hq) + (tateYtail a q hq - tateYtail b q hq)
        - (tateYterm (wOf q a) - tateYterm (wOf q b))
        - (tateYtail (wOf q a) q hq - tateYtail (wOf q b) q hq) := by ring
    rw [hkey]
    exact Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.add_mem _ dTfa dTga) dgw) dTgw
  rw [tateEps_eq a (wOf q a) q hq ha1, tateEps_eq b (wOf q b) q hq hb1]
  have hfinal : ((tateYtail a q hq - tateXterm (wOf q a) - tateXtail (wOf q a) q hq
          - tateYterm (wOf q a) - tateYtail (wOf q a) q hq + evalAdic (sigmaSeries 1) q hq)
        - a * (tateXtail a q hq + tateYtail a q hq - tateYterm (wOf q a)
          - tateYtail (wOf q a) q hq - evalAdic (sigmaSeries 1) q hq))
      - ((tateYtail b q hq - tateXterm (wOf q b) - tateXtail (wOf q b) q hq
          - tateYterm (wOf q b) - tateYtail (wOf q b) q hq + evalAdic (sigmaSeries 1) q hq)
        - b * (tateXtail b q hq + tateYtail b q hq - tateYterm (wOf q b)
          - tateYtail (wOf q b) q hq - evalAdic (sigmaSeries 1) q hq))
      = (tateYtail a q hq - tateYtail b q hq)
        - (tateXterm (wOf q a) - tateXterm (wOf q b))
        - (tateXtail (wOf q a) q hq - tateXtail (wOf q b) q hq)
        - (tateYterm (wOf q a) - tateYterm (wOf q b))
        - (tateYtail (wOf q a) q hq - tateYtail (wOf q b) q hq)
        - (a - b) * (tateXtail a q hq + tateYtail a q hq - tateYterm (wOf q a)
          - tateYtail (wOf q a) q hq - evalAdic (sigmaSeries 1) q hq)
        - b * ((tateXtail a q hq + tateYtail a q hq - tateYterm (wOf q a)
            - tateYtail (wOf q a) q hq - evalAdic (sigmaSeries 1) q hq)
          - (tateXtail b q hq + tateYtail b q hq - tateYterm (wOf q b)
            - tateYtail (wOf q b) q hq - evalAdic (sigmaSeries 1) q hq)) := by ring
  rw [hfinal]
  refine Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _
    (Ideal.sub_mem _ (Ideal.sub_mem _ dTga dfw) dTfw) dgw) dTgw) ?_) ?_
  · rw [pow_succ]
    exact Ideal.mul_mem_mul hab hE2a
  · exact Ideal.mul_mem_left _ _ dE2

/-! ## ★出典の紐付け(`.src`) -/

def wOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——単元 a の相方 q/a)",
    sectionId := "genell-def-3-3" }

def tateEps_diff_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——誤差項の差は 1 つ位が上がる)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
