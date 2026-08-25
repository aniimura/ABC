import ABC3.Found.GaloisRep.TateLinearQ

/-!
# Galois (G6) 第 271 ブロック —— **★★★★★★★★★環帯での全射性が閉じた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★到達点

> `x, y ∈ I` で `(x,y)` が `E_q` 上なら、`a, w ∈ I` が在って
> **`a·w = q`** かつ `X(a,w,q) = x`、`Y(a,w,q) = y`(`tate_surjective_annulus`)

★第 269 では制約 `a·w = q` を課さずに `(a,w)` を作った。本ブロックでそれを回復する。

## ★★★★★★★★`defect = 0` は `q` を決める

    defect(a,w,q) = 0  かつ  defect(a,w,a·w) = 0   ⟹   q = a·w

★逐次近似:`q − q' ∈ I^k` と仮定すると、曲線の式の差から

    a₆(q) − a₆(q') = (Y−Y')(Y+Y') + (X−X')Y + (Y−Y')X'
                     − (X−X')(X²+XX'+X'²) − (a₄−a₄')X − (X−X')a₄'

の右辺はすべて **`I^k` の元 × `I` の元** なので `I^{k+1}`。
一方 第 270 より `a₆(q) − a₆(q') + (q − q') ∈ I^{k+1}` なので **`q − q' ∈ I^{k+1}`**。
`IsHausdorff` で `q = q'`。

★★**`X, Y, a₄` がすべて `I` に入る**(環帯だから)ことが効いている。
`a₆` だけが `q` を 1 次で拾い、他は全部 1 つ位を上げる——それが「`q` が決まる」ことの中身。

## ★★`q` についての差

`X(a,w,q) − X(a,w,q')` は尾と `s₁(q)` の差だけ(`f(a)`, `f(w)` は `q` によらない)。
各項は `q^{n+1}u − q'^{n+1}u = (q^{n+1} − q'^{n+1})u` が `q − q'` で割れるので `I^k`。
`adicSum_diff_mem`(第 259)でそのまま和に上がる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateXpair_diff_q_mem`・`tateYpair_diff_q_mem` | ★★★★`q` についての差 |
| `tateXpair_mem`・`tateYpair_mem` | ★★環帯では座標も `I` の元 |
| `tate_q_unique` | ★★★★★★★★**`defect = 0` は `q` を決める** |
| `tate_surjective_annulus` | ★★★★★★★★★**環帯での全射性(制約つき)** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★`q` についての差 -/

theorem pow_succ_diff_mul_mem {k n : ℕ} {q q' u : R} (hqq : q - q' ∈ I ^ k) :
    q ^ (n + 1) * u - q' ^ (n + 1) * u ∈ I ^ k := by
  have h : q ^ (n + 1) * u - q' ^ (n + 1) * u = (q ^ (n + 1) - q' ^ (n + 1)) * u := by ring
  rw [h]
  exact Ideal.mul_mem_right _ _
    (Ideal.mem_of_dvd _ (sub_dvd_pow_sub_pow q q' (n + 1)) hqq)

theorem tateXtail_diff_q_mem [IsAdicComplete I R] {k : ℕ} (u q q' : R) (hq : q ∈ I)
    (hq' : q' ∈ I) (hqq : q - q' ∈ I ^ k) :
    tateXtail u q hq - tateXtail u q' hq' ∈ I ^ k := by
  refine adicSum_diff_mem _ _ _ _ k fun n => ?_
  exact tateXterm_diff_mem' (pow_succ_mul_mem' hq n u) (pow_succ_mul_mem' hq' n u)
    (pow_succ_diff_mul_mem hqq)

theorem tateYtail_diff_q_mem [IsAdicComplete I R] {k : ℕ} (u q q' : R) (hq : q ∈ I)
    (hq' : q' ∈ I) (hqq : q - q' ∈ I ^ k) :
    tateYtail u q hq - tateYtail u q' hq' ∈ I ^ k := by
  refine adicSum_diff_mem _ _ _ _ k fun n => ?_
  exact tateYterm_diff_mem' (pow_succ_mul_mem' hq n u) (pow_succ_mul_mem' hq' n u)
    (pow_succ_diff_mul_mem hqq)

theorem tateXpair_diff_q_mem [IsAdicComplete I R] {k : ℕ} (a w q q' : R) (hq : q ∈ I)
    (hq' : q' ∈ I) (hqq : q - q' ∈ I ^ k) :
    tateXpair a w q hq - tateXpair a w q' hq' ∈ I ^ k := by
  have hkey : tateXpair a w q hq - tateXpair a w q' hq'
      = (tateXtail a q hq - tateXtail a q' hq') + (tateXtail w q hq - tateXtail w q' hq')
        + (-2) * (evalAdic (sigmaSeries 1) q hq - evalAdic (sigmaSeries 1) q' hq') := by
    rw [tateXpair, tateXpair]; ring
  rw [hkey]
  exact Ideal.add_mem _ (Ideal.add_mem _ (tateXtail_diff_q_mem a q q' hq hq' hqq)
    (tateXtail_diff_q_mem w q q' hq hq' hqq))
    (Ideal.mul_mem_left _ _ (evalAdic_sub_mem (sigmaSeries 1) hq hq' k hqq))

theorem tateYpair_diff_q_mem [IsAdicComplete I R] {k : ℕ} (a w q q' : R) (hq : q ∈ I)
    (hq' : q' ∈ I) (hqq : q - q' ∈ I ^ k) :
    tateYpair a w q hq - tateYpair a w q' hq' ∈ I ^ k := by
  have hkey : tateYpair a w q hq - tateYpair a w q' hq'
      = (tateYtail a q hq - tateYtail a q' hq') + (-1) * (tateXtail w q hq - tateXtail w q' hq')
        + (-1) * (tateYtail w q hq - tateYtail w q' hq')
        + (evalAdic (sigmaSeries 1) q hq - evalAdic (sigmaSeries 1) q' hq') := by
    rw [tateYpair, tateYpair]; ring
  rw [hkey]
  exact Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _
    (tateYtail_diff_q_mem a q q' hq hq' hqq)
    (Ideal.mul_mem_left _ _ (tateXtail_diff_q_mem w q q' hq hq' hqq)))
    (Ideal.mul_mem_left _ _ (tateYtail_diff_q_mem w q q' hq hq' hqq)))
    (evalAdic_sub_mem (sigmaSeries 1) hq hq' k hqq)

/-! ## ★★環帯では座標も `I` の元 -/

theorem tateXpair_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (ha : a ∈ I) (hw : w ∈ I) :
    tateXpair a w q hq ∈ I := by
  have h := tateXpair_sub_mem a w q hq hw
  have h2 := Ideal.add_mem I h (tateXterm_mem_one ha)
  simpa using h2

theorem tateYpair_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (ha : a ∈ I) (hw : w ∈ I) :
    tateYpair a w q hq ∈ I := by
  have h := tateYpair_sub_tateYterm_mem a w q hq hw
  have h2 := Ideal.add_mem I h (tateYterm_mem_one ha)
  simpa using h2

/-! ## ★★★★★★★★`defect = 0` は `q` を決める -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**`defect = 0` は `q` を決める**——`a₆` が `q` を 1 次で拾うから。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_q_unique [IsAdicComplete I R] (a w q q' : R) (hq : q ∈ I) (hq' : q' ∈ I)
    (ha : a ∈ I) (hw : w ∈ I)
    (hd : tateYpair a w q hq ^ 2 + tateXpair a w q hq * tateYpair a w q hq
      = tateXpair a w q hq ^ 3 + (tateCurveAt q hq).a₄ * tateXpair a w q hq
        + (tateCurveAt q hq).a₆)
    (hd' : tateYpair a w q' hq' ^ 2 + tateXpair a w q' hq' * tateYpair a w q' hq'
      = tateXpair a w q' hq' ^ 3 + (tateCurveAt q' hq').a₄ * tateXpair a w q' hq'
        + (tateCurveAt q' hq').a₆) : q = q' := by
  have hX := tateXpair_mem a w q hq ha hw
  have hX' := tateXpair_mem a w q' hq' ha hw
  have hY := tateYpair_mem a w q hq ha hw
  have hY' := tateYpair_mem a w q' hq' ha hw
  have h4' := tateCurveAt_a4_mem (I := I) q' hq'
  have hbr : tateXpair a w q hq ^ 2 + tateXpair a w q hq * tateXpair a w q' hq'
      + tateXpair a w q' hq' ^ 2 ∈ I := by
    refine Ideal.add_mem _ (Ideal.add_mem _ ?_ ?_) ?_
    · rw [sq]; exact Ideal.mul_mem_right _ _ hX
    · exact Ideal.mul_mem_right _ _ hX
    · rw [sq]; exact Ideal.mul_mem_right _ _ hX'
  have hstep : ∀ k : ℕ, q - q' ∈ I ^ k := by
    intro k
    induction k with
    | zero => simp
    | succ n ih =>
      have dX := tateXpair_diff_q_mem a w q q' hq hq' ih
      have dY := tateYpair_diff_q_mem a w q q' hq hq' ih
      have d4 := tateCurveAt_a4_sub_mem hq hq' n ih
      have hkey : (tateCurveAt q hq).a₆ - (tateCurveAt q' hq').a₆
          = (tateYpair a w q hq - tateYpair a w q' hq')
              * (tateYpair a w q hq + tateYpair a w q' hq')
            + (tateXpair a w q hq - tateXpair a w q' hq') * tateYpair a w q hq
            + (tateYpair a w q hq - tateYpair a w q' hq') * tateXpair a w q' hq'
            - (tateXpair a w q hq - tateXpair a w q' hq')
              * (tateXpair a w q hq ^ 2
                + tateXpair a w q hq * tateXpair a w q' hq'
                + tateXpair a w q' hq' ^ 2)
            - ((tateCurveAt q hq).a₄ - (tateCurveAt q' hq').a₄) * tateXpair a w q hq
            - (tateXpair a w q hq - tateXpair a w q' hq') * (tateCurveAt q' hq').a₄ := by
        linear_combination hd' - hd
      have ha6 : (tateCurveAt q hq).a₆ - (tateCurveAt q' hq').a₆ ∈ I ^ (n + 1) := by
        rw [hkey, pow_succ]
        exact Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.add_mem _
          (Ideal.add_mem _ (Ideal.mul_mem_mul dY (Ideal.add_mem _ hY hY'))
          (Ideal.mul_mem_mul dX hY)) (Ideal.mul_mem_mul dY hX'))
          (Ideal.mul_mem_mul dX hbr)) (Ideal.mul_mem_mul d4 hX)) (Ideal.mul_mem_mul dX h4')
      have ha6' := tateCurveAt_a6_sub_linear hq hq' n ih
      have hfin : q - q' = ((tateCurveAt q hq).a₆ - (tateCurveAt q' hq').a₆ + (q - q'))
          - ((tateCurveAt q hq).a₆ - (tateCurveAt q' hq').a₆) := by ring
      rw [hfin]
      exact Ideal.sub_mem _ ha6' ha6
  exact sub_eq_zero.1 (eq_zero_of_mem_pow (I := I) hstep)

/-! ## ★★★★★★★★★環帯での全射性 -/

set_option maxHeartbeats 1200000 in
/-- ★★★★★★★★★**環帯での全射性(制約つき)**——葉 (e) の環帯の領域。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_surjective_annulus [IsAdicComplete I R] (x y q : R) (hq : q ∈ I)
    (hx : x ∈ I) (hy : y ∈ I)
    (he : y ^ 2 + x * y = x ^ 3 + (tateCurveAt q hq).a₄ * x + (tateCurveAt q hq).a₆) :
    ∃ a ∈ I, ∃ w ∈ I, a * w = q ∧ tateXpair a w q hq = x ∧ tateYpair a w q hq = y := by
  obtain ⟨a, ha, w, hw, hX, hY⟩ := exists_tate_annulus x y q hq hx hy
  have haw : a * w ∈ I := Ideal.mul_mem_right _ _ ha
  have hd : tateYpair a w q hq ^ 2 + tateXpair a w q hq * tateYpair a w q hq
      = tateXpair a w q hq ^ 3 + (tateCurveAt q hq).a₄ * tateXpair a w q hq
        + (tateCurveAt q hq).a₆ := by
    rw [hX, hY]; exact he
  have hd' := tate_equation a w (a * w) haw rfl (isUnit_one_sub ha) (isUnit_one_sub hw)
  have hqeq : q = a * w := tate_q_unique a w q (a * w) hq haw ha hw hd hd'
  exact ⟨a, ha, w, hw, hqeq.symm, hX, hY⟩

/-! ## ★出典の紐付け(`.src`) -/

def tate_q_unique.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——defect = 0 は q を決める)",
    sectionId := "genell-def-3-3" }

def tate_surjective_annulus.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——環帯での全射性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
