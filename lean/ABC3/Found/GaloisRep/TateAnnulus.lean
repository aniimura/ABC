import ABC3.Found.GaloisRep.AdicContraction2

/-!
# Galois (G6) 第 269 ブロック —— **★★★★★★★★環帯での座標の全射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★到達点

> `x, y ∈ I` なら `a, w ∈ I` が在って `X(a,w,q) = x`、`Y(a,w,q) = y`
> (`exists_tate_annulus`)

★**制約 `a·w = q` はまだ課していない**——それは次のブロックで
「`defect` の `q` についての主要項が `q − aw`」から回復させる(§9-581)。

## ★★★★★2 つの主要項

環帯(`a, w ∈ I`)では

    X + Y ≡ a,     Y ≡ −w        (誤差はどちらも `I` の元)

だから、誤差を `δ_S := (X+Y) − a`、`δ_Y := Y + w` と置いて

    a = (x + y) − δ_S(a,w),      w = −y + δ_Y(a,w)

の 2 変数不動点(第 268)を取る。不動点では `X + Y = x + y` と `Y = y`、
したがって `X = x` が出る。

## ★★★誤差が 1 つ位を上げること

`δ_S = (f(a)+g(a)−a) + Tf(a) + Tg(a) − g(w) − Tg(w) − s₁` で、`s₁(q)` は
両側で共通なので**差では消える**。残りは

| 項 | 差の位 |
|---|---|
| `f(a)+g(a)−a = a(inv(1−a)³ − 1)` | `I^{k+1}`(★下記) |
| 尾 `Tf(a), Tg(a)` | `I^{k+1}`(`q ∈ I` を 1 つ含む) |
| `g(w)`, `Tg(w)` | `I^{k+1}` |

★**`f + g` の主要項がちょうど `a` である**(第 262 の `mul_tateXYterm`)ことが効いて、
`f(a)+g(a)−a = a·(inv(1−a)³ − 1)` と**積の形**に書ける。差は
`(a−a')·h(a) + (h(a)−h(a'))·a'` に分けるだけで `I^{k+1}` に入る。

`δ_Y` も同様で、`f(w) − w` の差が 1 つ位を上げること(第 258 の `tateXterm_diff_mem`)を使う。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `inverse_pow_sub_one_mem` | ★★★`inv((1−a)ⁿ) − 1 ∈ I` |
| `tateXYterm_sub_self_eq` | ★★★★`f(a)+g(a)−a = a(inv(1−a)³−1)` |
| `tateDeltaS`・`tateDeltaY` | ★★★環帯での 2 つの誤差 |
| `tateDeltaS_diff_mem`・`tateDeltaY_diff_mem` | ★★★★★誤差は 1 つ位を上げる |
| `exists_tate_annulus` | ★★★★★★★★**環帯での座標の全射性** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★逆元の主要項 -/

theorem one_sub_pow_sub_one_mem {a : R} (ha : a ∈ I) (n : ℕ) : (1 - a) ^ n - 1 ∈ I := by
  refine Ideal.Quotient.eq.1 ?_
  have h : (Ideal.Quotient.mk I) a = 0 := Ideal.Quotient.eq_zero_iff_mem.2 ha
  simp [map_pow, map_sub, map_one, h]

theorem inverse_pow_sub_one_mem [IsAdicComplete I R] {a : R} (ha : a ∈ I) (n : ℕ) :
    Ring.inverse ((1 - a) ^ n) - 1 ∈ I := by
  have hu : IsUnit ((1 - a) ^ n) := (isUnit_one_sub ha).pow n
  have h1 : Ring.inverse ((1 - a) ^ n) * (1 - a) ^ n = 1 := Ring.inverse_mul_cancel _ hu
  have hkey : Ring.inverse ((1 - a) ^ n) - 1
      = Ring.inverse ((1 - a) ^ n) * (1 - (1 - a) ^ n) := by
    rw [mul_sub, mul_one, h1]
  rw [hkey]
  refine Ideal.mul_mem_left _ _ ?_
  have h2 := neg_mem (one_sub_pow_sub_one_mem ha n)
  rwa [neg_sub] at h2

/-! ## ★★★★`f + g` の主要項 -/

/-- ★★★★**`f(a) + g(a) − a = a·(inv(1−a)³ − 1)`**——積の形に書ける。 -/
theorem tateXYterm_sub_self_eq [IsAdicComplete I R] {a : R} (ha : a ∈ I) :
    tateXterm a + tateYterm a - a = a * (Ring.inverse ((1 - a) ^ 3) - 1) := by
  have hu : IsUnit (1 - a) := isUnit_one_sub ha
  have h1 : Ring.inverse ((1 - a) ^ 3) * (1 - a) ^ 3 = 1 :=
    Ring.inverse_mul_cancel _ (hu.pow 3)
  have hsum := mul_tateXYterm hu
  have hkey : (1 - a) ^ 3 * ((tateXterm a + tateYterm a - a)
      - a * (Ring.inverse ((1 - a) ^ 3) - 1)) = 0 := by
    have h2 : (1 - a) ^ 3 * (a * (Ring.inverse ((1 - a) ^ 3) - 1))
        = a * ((1 - a) ^ 3 * Ring.inverse ((1 - a) ^ 3)) - a * (1 - a) ^ 3 := by ring
    rw [mul_sub, h2, Ring.mul_inverse_cancel _ (hu.pow 3), mul_one, mul_sub, hsum]
    ring
  have h3 := congrArg (fun z => Ring.inverse ((1 - a) ^ 3) * z) hkey
  simp only [← mul_assoc, h1, one_mul, mul_zero] at h3
  exact sub_eq_zero.1 h3

theorem tateXYterm_sub_self_mem [IsAdicComplete I R] {a : R} (ha : a ∈ I) :
    tateXterm a + tateYterm a - a ∈ I := by
  rw [tateXYterm_sub_self_eq ha]
  exact Ideal.mul_mem_left _ _ (inverse_pow_sub_one_mem ha 3)

theorem tateXYterm_sub_self_diff_mem [IsAdicComplete I R] {k : ℕ} {a b : R}
    (ha : a ∈ I) (hb : b ∈ I) (hab : a - b ∈ I ^ k) :
    (tateXterm a + tateYterm a - a) - (tateXterm b + tateYterm b - b) ∈ I ^ (k + 1) := by
  rw [tateXYterm_sub_self_eq ha, tateXYterm_sub_self_eq hb]
  have hcube : (1 - a) ^ 3 - (1 - b) ^ 3 ∈ I ^ k := by
    have h : (1 - a) ^ 3 - (1 - b) ^ 3
        = (b - a) * ((1 - a) ^ 2 + (1 - a) * (1 - b) + (1 - b) ^ 2) := by ring
    rw [h]
    refine Ideal.mul_mem_right _ _ ?_
    have h2 := neg_mem hab
    rwa [neg_sub] at h2
  have hinv : Ring.inverse ((1 - a) ^ 3) - Ring.inverse ((1 - b) ^ 3) ∈ I ^ k :=
    inverse_diff_mem ((isUnit_one_sub ha).pow 3) ((isUnit_one_sub hb).pow 3) hcube
  have hkey : a * (Ring.inverse ((1 - a) ^ 3) - 1) - b * (Ring.inverse ((1 - b) ^ 3) - 1)
      = (a - b) * (Ring.inverse ((1 - a) ^ 3) - 1)
        + (Ring.inverse ((1 - a) ^ 3) - Ring.inverse ((1 - b) ^ 3)) * b := by ring
  rw [hkey, pow_succ]
  exact Ideal.add_mem _ (Ideal.mul_mem_mul hab (inverse_pow_sub_one_mem ha 3))
    (Ideal.mul_mem_mul hinv hb)

/-! ## ★★★環帯での 2 つの誤差 -/

/-- ★★★環帯での誤差 `δ_S = (X + Y) − a`。 -/
noncomputable def tateDeltaS [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) : R :=
  tateXpair a w q hq + tateYpair a w q hq - a

/-- ★★★環帯での誤差 `δ_Y = Y + w`。 -/
noncomputable def tateDeltaY [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) : R :=
  tateYpair a w q hq + w

theorem tateDeltaS_eq [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) :
    tateDeltaS a w q hq = (tateXterm a + tateYterm a - a) + tateXtail a q hq + tateYtail a q hq
      - tateYterm w - tateYtail w q hq - evalAdic (sigmaSeries 1) q hq := by
  rw [tateDeltaS, tateXpair, tateYpair]; ring

theorem tateDeltaY_eq [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) :
    tateDeltaY a w q hq = tateYterm a + tateYtail a q hq - (tateXterm w - w)
      - tateXtail w q hq - tateYterm w - tateYtail w q hq + evalAdic (sigmaSeries 1) q hq := by
  rw [tateDeltaY, tateYpair]; ring

theorem tateXterm_sub_self_mem_one [IsAdicComplete I R] {w : R} (hw : w ∈ I) :
    tateXterm w - w ∈ I := by
  have h := tateXterm_sub_self_mem (I := I) (k := 1) (by rwa [pow_one]) hw
  have h2 : (2 : ℕ) * 1 = 2 := by norm_num
  rw [h2] at h
  have h3 : I ^ 2 ≤ I ^ 1 := Ideal.pow_le_pow_right (by omega)
  rw [pow_one] at h3
  exact h3 h

theorem tateDeltaS_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (ha : a ∈ I) (hw : w ∈ I) :
    tateDeltaS a w q hq ∈ I := by
  rw [tateDeltaS_eq]
  exact Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.add_mem _ (Ideal.add_mem _
    (tateXYterm_sub_self_mem ha) (tateXtail_mem a q hq)) (tateYtail_mem a q hq))
    (tateYterm_mem_one hw)) (tateYtail_mem w q hq))
    (evalAdic_mem (sigmaSeries 1) (by simp [sigmaSeries]) q hq)

theorem tateDeltaY_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (ha : a ∈ I) (hw : w ∈ I) :
    tateDeltaY a w q hq ∈ I := by
  rw [tateDeltaY_eq]
  exact Ideal.add_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _
    (Ideal.add_mem _ (tateYterm_mem_one ha) (tateYtail_mem a q hq))
    (tateXterm_sub_self_mem_one hw)) (tateXtail_mem w q hq))
    (tateYterm_mem_one hw)) (tateYtail_mem w q hq))
    (evalAdic_mem (sigmaSeries 1) (by simp [sigmaSeries]) q hq)

set_option maxHeartbeats 1200000 in
theorem tateDeltaS_diff_mem [IsAdicComplete I R] {k : ℕ} (a w a' w' q : R) (hq : q ∈ I)
    (ha : a ∈ I) (hw : w ∈ I) (ha' : a' ∈ I) (hw' : w' ∈ I)
    (hda : a - a' ∈ I ^ k) (hdw : w - w' ∈ I ^ k) :
    tateDeltaS a w q hq - tateDeltaS a' w' q hq ∈ I ^ (k + 1) := by
  have hq1 : q ∈ I ^ 1 := by rwa [pow_one]
  have d1 := tateXYterm_sub_self_diff_mem ha ha' hda
  have d2 := tateXtail_diff_mem (m := 1) a a' q hq hq1 hda
  have d3 := tateYtail_diff_mem (m := 1) a a' q hq hq1 hda
  have d4 := tateYterm_diff_mem hw hw' hdw
  have d5 := tateYtail_diff_mem (m := 1) w w' q hq hq1 hdw
  rw [tateDeltaS_eq, tateDeltaS_eq]
  have hkey : ((tateXterm a + tateYterm a - a) + tateXtail a q hq + tateYtail a q hq
        - tateYterm w - tateYtail w q hq - evalAdic (sigmaSeries 1) q hq)
      - ((tateXterm a' + tateYterm a' - a') + tateXtail a' q hq + tateYtail a' q hq
        - tateYterm w' - tateYtail w' q hq - evalAdic (sigmaSeries 1) q hq)
      = ((tateXterm a + tateYterm a - a) - (tateXterm a' + tateYterm a' - a'))
        + (tateXtail a q hq - tateXtail a' q hq) + (tateYtail a q hq - tateYtail a' q hq)
        - (tateYterm w - tateYterm w') - (tateYtail w q hq - tateYtail w' q hq) := by ring
  rw [hkey]
  exact Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.add_mem _ (Ideal.add_mem _ d1 d2) d3) d4) d5

set_option maxHeartbeats 1200000 in
theorem tateDeltaY_diff_mem [IsAdicComplete I R] {k : ℕ} (a w a' w' q : R) (hq : q ∈ I)
    (ha : a ∈ I) (hw : w ∈ I) (ha' : a' ∈ I) (hw' : w' ∈ I)
    (hda : a - a' ∈ I ^ k) (hdw : w - w' ∈ I ^ k) :
    tateDeltaY a w q hq - tateDeltaY a' w' q hq ∈ I ^ (k + 1) := by
  have hq1 : q ∈ I ^ 1 := by rwa [pow_one]
  have e1 := tateYterm_diff_mem ha ha' hda
  have e2 := tateYtail_diff_mem (m := 1) a a' q hq hq1 hda
  have e3 := tateXterm_diff_mem hw hw' hdw
  have e4 := tateXtail_diff_mem (m := 1) w w' q hq hq1 hdw
  have e5 := tateYterm_diff_mem hw hw' hdw
  have e6 := tateYtail_diff_mem (m := 1) w w' q hq hq1 hdw
  rw [tateDeltaY_eq, tateDeltaY_eq]
  have hkey : (tateYterm a + tateYtail a q hq - (tateXterm w - w) - tateXtail w q hq
        - tateYterm w - tateYtail w q hq + evalAdic (sigmaSeries 1) q hq)
      - (tateYterm a' + tateYtail a' q hq - (tateXterm w' - w') - tateXtail w' q hq
        - tateYterm w' - tateYtail w' q hq + evalAdic (sigmaSeries 1) q hq)
      = (tateYterm a - tateYterm a') + (tateYtail a q hq - tateYtail a' q hq)
        - (tateXterm w - tateXterm w' - (w - w'))
        - (tateXtail w q hq - tateXtail w' q hq)
        - (tateYterm w - tateYterm w') - (tateYtail w q hq - tateYtail w' q hq) := by ring
  rw [hkey]
  exact Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.sub_mem _
    (Ideal.add_mem _ e1 e2) e3) e4) e5) e6

/-! ## ★★★★★★★★環帯での全射性 -/

set_option maxHeartbeats 1200000 in
/-- ★★★★★★★★**環帯での座標の全射性**——制約 `a·w = q` はまだ課さない。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_tate_annulus [IsAdicComplete I R] (x y q : R) (hq : q ∈ I)
    (hx : x ∈ I) (hy : y ∈ I) :
    ∃ a ∈ I, ∃ w ∈ I, tateXpair a w q hq = x ∧ tateYpair a w q hq = y := by
  have hFI : ∀ a ∈ I, ∀ w ∈ I, (x + y) - tateDeltaS a w q hq ∈ I := fun a ha w hw =>
    Ideal.sub_mem _ (Ideal.add_mem _ hx hy) (tateDeltaS_mem a w q hq ha hw)
  have hGI : ∀ a ∈ I, ∀ w ∈ I, -y + tateDeltaY a w q hq ∈ I := fun a ha w hw =>
    Ideal.add_mem _ (neg_mem hy) (tateDeltaY_mem a w q hq ha hw)
  have hconF : ∀ a ∈ I, ∀ w ∈ I, ∀ a' ∈ I, ∀ w' ∈ I, ∀ k : ℕ,
      a - a' ∈ I ^ k → w - w' ∈ I ^ k →
      ((x + y) - tateDeltaS a w q hq) - ((x + y) - tateDeltaS a' w' q hq) ∈ I ^ (k + 1) := by
    intro a ha w hw a' ha' w' hw' k hda hdw
    have h : ((x + y) - tateDeltaS a w q hq) - ((x + y) - tateDeltaS a' w' q hq)
        = -(tateDeltaS a w q hq - tateDeltaS a' w' q hq) := by ring
    rw [h]
    exact neg_mem (tateDeltaS_diff_mem a w a' w' q hq ha hw ha' hw' hda hdw)
  have hconG : ∀ a ∈ I, ∀ w ∈ I, ∀ a' ∈ I, ∀ w' ∈ I, ∀ k : ℕ,
      a - a' ∈ I ^ k → w - w' ∈ I ^ k →
      (-y + tateDeltaY a w q hq) - (-y + tateDeltaY a' w' q hq) ∈ I ^ (k + 1) := by
    intro a ha w hw a' ha' w' hw' k hda hdw
    have h : (-y + tateDeltaY a w q hq) - (-y + tateDeltaY a' w' q hq)
        = tateDeltaY a w q hq - tateDeltaY a' w' q hq := by ring
    rw [h]
    exact tateDeltaY_diff_mem a w a' w' q hq ha hw ha' hw' hda hdw
  obtain ⟨a, ha, w, hw, hF, hG⟩ := exists_fixedPoint_of_contraction₂
    (fun a w => (x + y) - tateDeltaS a w q hq) (fun a w => -y + tateDeltaY a w q hq)
    hFI hGI hconF hconG
  have hY : tateYpair a w q hq = y := by
    have h : tateDeltaY a w q hq = w + y := by linear_combination hG
    rw [tateDeltaY] at h
    linear_combination h
  have hX : tateXpair a w q hq = x := by
    have h : tateDeltaS a w q hq = (x + y) - a := by linear_combination -hF
    rw [tateDeltaS] at h
    linear_combination h - hY
  exact ⟨a, ha, w, hw, hX, hY⟩

/-! ## ★出典の紐付け(`.src`) -/

def tateXYterm_sub_self_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——f+g の主要項は a)",
    sectionId := "genell-def-3-3" }

def exists_tate_annulus.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——環帯での座標の全射性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
