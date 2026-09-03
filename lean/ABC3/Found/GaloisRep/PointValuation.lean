/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfBaseChange
import ABC3.Meta.Claim

/-!
# 第 1070 ブロック —— **点の座標の付値**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か

第 1069 の測定により、`Lemma 3.5` に残るのは
`isIntegral_veluQuotientFull_of_coprime` の**良い素点側**ただ 1 つで、
その中身は「位数が剰余標数と素な捻れ点は極小モデルで整座標をもつ」である。

☆古典的な道は**還元の核の濾過**（形式群）であり、その**第 1 段**が本ファイルである:

> `W` が `p` で整、`(x, y)` がその上の点、`v_p(x) < 0` ならば
> `2·v_p(y) = 3·v_p(x)`。すなわち `v_p(x) = −2m`・`v_p(y) = −3m`。

★これは Weierstrass 方程式の**両辺の主要項が `y²` と `x³` しかない**ことから出る、
純粋な超距離不等式の帰結である。曲線論も還元も要らない。

## ★★★★`ValAtLeast` —— 零を許す `valAdd ≥ n`

`valAdd p` は `Lˣ` の上でしか定義されないので、`a₂ x² + a₄ x + a₆` のように
**個々の項が零になりうる和**を扱えない。そこで

    ValAtLeast p n z  :⟺  ∀ (hz : z ≠ 0), n ≤ valAdd p ⟨z⟩

を置く。`z = 0` では**空虚に真**であり、和・積・単元倍で閉じる。
☆これで超距離不等式の帳簿が `omega` に載る形になる。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve ABC3.Meta

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★`ValAtLeast` の API -/

/-- ★★`valAdd p z ≥ n`（`z = 0` も許す形）。 -/
def ValAtLeast (p : HeightOneSpectrum (𝓞 L)) (n : ℤ) (z : L) : Prop :=
  ∀ hz : z ≠ 0, n ≤ valAdd p (Units.mk0 z hz)

theorem valAtLeast_zero (p : HeightOneSpectrum (𝓞 L)) (n : ℤ) : ValAtLeast p n 0 :=
  fun hz => absurd rfl hz

theorem valAtLeast_mono {p : HeightOneSpectrum (𝓞 L)} {m n : ℤ} {z : L}
    (h : m ≤ n) (hz : ValAtLeast p n z) : ValAtLeast p m z :=
  fun h0 => le_trans h (hz h0)

/-- ★★単元の `valAdd` はそのまま `ValAtLeast` になる。 -/
theorem valAtLeast_unit (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) :
    ValAtLeast p (valAdd p x) (x : L) := fun hz =>
  le_of_eq (valAdd_eq_of_valuation_eq p x (Units.mk0 (x : L) hz) rfl)

/-- ★★付値環の元は `ValAtLeast p 0`。 -/
theorem valAtLeast_of_mem {p : HeightOneSpectrum (𝓞 L)} {z : L}
    (h : z ∈ primeSubring p) : ValAtLeast p 0 z := by
  intro hz
  rw [valAdd_nonneg_iff]
  exact (mem_primeSubring_iff p z).1 h

theorem valAtLeast_neg {p : HeightOneSpectrum (𝓞 L)} {n : ℤ} {z : L}
    (h : ValAtLeast p n z) : ValAtLeast p n (-z) := by
  intro hz
  have hz0 : z ≠ 0 := fun hc => hz (by simp [hc])
  refine le_trans (h hz0) (le_of_eq (valAdd_eq_of_valuation_eq p _ _ ?_))
  simp

/-- ★★★★**超距離不等式**——和は最小値を下回らない。 -/
theorem valAtLeast_add {p : HeightOneSpectrum (𝓞 L)} {n : ℤ} {a b : L}
    (ha : ValAtLeast p n a) (hb : ValAtLeast p n b) : ValAtLeast p n (a + b) := by
  intro hab
  rcases eq_or_ne a 0 with rfl | ha0
  · refine le_trans (hb (by simpa using hab)) (le_of_eq (valAdd_eq_of_valuation_eq p _ _ ?_))
    simp
  rcases eq_or_ne b 0 with rfl | hb0
  · refine le_trans (ha (by simpa using hab)) (le_of_eq (valAdd_eq_of_valuation_eq p _ _ ?_))
    simp
  have h1 := ha ha0
  have h2 := hb hb0
  have hv : (p.valuation L) (a + b) ≤ max ((p.valuation L) a) ((p.valuation L) b) :=
    Valuation.map_add _ _ _
  rcases le_total (valAdd p (Units.mk0 a ha0)) (valAdd p (Units.mk0 b hb0)) with hle | hle
  · have hba : (p.valuation L) b ≤ (p.valuation L) a := (valAdd_le_iff p _ _).1 hle
    have hstep : (p.valuation L) (a + b) ≤ (p.valuation L) a := by
      refine le_trans hv ?_
      rw [max_eq_left hba]
    exact le_trans h1 ((valAdd_le_iff p (Units.mk0 a ha0) (Units.mk0 (a + b) hab)).2 hstep)
  · have hab' : (p.valuation L) a ≤ (p.valuation L) b := (valAdd_le_iff p _ _).1 hle
    have hstep : (p.valuation L) (a + b) ≤ (p.valuation L) b := by
      refine le_trans hv ?_
      rw [max_eq_right hab']
    exact le_trans h2 ((valAdd_le_iff p (Units.mk0 b hb0) (Units.mk0 (a + b) hab)).2 hstep)

/-- ★★★積の `ValAtLeast` は足し算になる。 -/
theorem valAtLeast_mul {p : HeightOneSpectrum (𝓞 L)} {m n : ℤ} {a b : L}
    (ha : ValAtLeast p m a) (hb : ValAtLeast p n b) : ValAtLeast p (m + n) (a * b) := by
  intro hab
  have ha0 : a ≠ 0 := fun h => hab (by simp [h])
  have hb0 : b ≠ 0 := fun h => hab (by simp [h])
  have heq : valAdd p (Units.mk0 (a * b) hab)
      = valAdd p (Units.mk0 a ha0) + valAdd p (Units.mk0 b hb0) := by
    rw [← valAdd_mul p (Units.mk0 a ha0) (Units.mk0 b hb0)]
    exact valAdd_eq_of_valuation_eq p _ _ (by simp)
  rw [heq]
  exact add_le_add (ha ha0) (hb hb0)

/-- ★★★★**厳密に小さい項があれば和はゼロにならない**。 -/
theorem add_ne_zero_of_valAdd_lt {p : HeightOneSpectrum (𝓞 L)} {n : ℤ} {a b : L}
    (ha : a ≠ 0) (hb : ValAtLeast p n b) (hlt : valAdd p (Units.mk0 a ha) < n) :
    a + b ≠ 0 := by
  intro hc
  have hb' : b = -a := by linear_combination hc
  subst hb'
  have hbne : (-a) ≠ 0 := by simpa using ha
  have h2 := hb hbne
  have h3 : valAdd p (Units.mk0 (-a) hbne) = valAdd p (Units.mk0 a ha) :=
    valAdd_eq_of_valuation_eq p _ _ (by simp)
  omega

/-- ★★★★**厳密に小さい項が和の `valAdd` を決める**。 -/
theorem valAdd_add_eq_of_lt {p : HeightOneSpectrum (𝓞 L)} {n : ℤ} {a b : L}
    (ha : a ≠ 0) (hab : a + b ≠ 0) (hb : ValAtLeast p n b)
    (hlt : valAdd p (Units.mk0 a ha) < n) :
    valAdd p (Units.mk0 (a + b) hab) = valAdd p (Units.mk0 a ha) := by
  rcases eq_or_ne b 0 with rfl | hb0
  · exact valAdd_eq_of_valuation_eq p _ _ (by simp)
  have h2 := hb hb0
  have hlt' : valAdd p (Units.mk0 a ha) < valAdd p (Units.mk0 b hb0) := lt_of_lt_of_le hlt h2
  have hle : (p.valuation L) b ≤ (p.valuation L) a := (valAdd_le_iff p _ _).1 hlt'.le
  have hne : (p.valuation L) b ≠ (p.valuation L) a := by
    intro h
    have : valAdd p (Units.mk0 b hb0) ≤ valAdd p (Units.mk0 a ha) :=
      (valAdd_le_iff p (Units.mk0 b hb0) (Units.mk0 a ha)).2 (le_of_eq h.symm)
    omega
  have hvlt : (p.valuation L) b < (p.valuation L) a := lt_of_le_of_ne hle hne
  exact valAdd_eq_of_valuation_eq p _ _
    (by simpa using Valuation.map_add_eq_of_lt_left (p.valuation L) hvlt)

/-- ☆**値が等しければ `valAdd` も等しい**——`Units.mk0` の中の `rw` を避ける道具。 -/
theorem valAdd_congr (p : HeightOneSpectrum (𝓞 L)) {a b : L} (ha : a ≠ 0) (hb : b ≠ 0)
    (h : a = b) : valAdd p (Units.mk0 a ha) = valAdd p (Units.mk0 b hb) :=
  valAdd_eq_of_valuation_eq p _ _ (by show (p.valuation L) a = (p.valuation L) b; rw [h])

/-! ## ★★★★★★★★★★★★★★★★第 1070 —— `2 v(y) = 3 v(x)` -/

section Point

variable (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
  [WeierstrassCurve.IsIntegral (primeSubring p) W]

/-- ☆整モデルの係数はすべて `ValAtLeast p 0`。 -/
private theorem valAtLeast_coeffs :
    ValAtLeast p 0 W.a₁ ∧ ValAtLeast p 0 W.a₂ ∧ ValAtLeast p 0 W.a₃ ∧
      ValAtLeast p 0 W.a₄ ∧ ValAtLeast p 0 W.a₆ := by
  obtain ⟨h1, h2, h3, h4, h6⟩ := mem_primeSubring_of_isIntegral p W
  exact ⟨valAtLeast_of_mem h1, valAtLeast_of_mem h2, valAtLeast_of_mem h3,
    valAtLeast_of_mem h4, valAtLeast_of_mem h6⟩

/-- ★★★★★★★★**方程式の右辺の `valAdd` は `3 v(x)`**（`v(x) < 0` のとき）。

☆`x³` が他のすべての項より厳密に深いので、和の `valAdd` は `x³` のそれで決まる。 -/
theorem valAdd_rhs_eq {x : L} (hx : x ≠ 0) (hneg : valAdd p (Units.mk0 x hx) < 0) :
    ∃ hne : x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆ ≠ 0,
      valAdd p (Units.mk0 (x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆) hne)
        = 3 * valAdd p (Units.mk0 x hx) := by
  obtain ⟨_, h2, _, h4, h6⟩ := valAtLeast_coeffs p W
  set a : ℤ := valAdd p (Units.mk0 x hx) with ha_def
  have hxA : ValAtLeast p a x := valAtLeast_unit p (Units.mk0 x hx)
  have hx3ne : x ^ 3 ≠ 0 := pow_ne_zero _ hx
  have hx3 : valAdd p (Units.mk0 (x ^ 3) hx3ne) = 3 * a := by
    have hmul : valAdd p (Units.mk0 x hx ^ 3) = ((3 : ℕ) : ℤ) * a := valAdd_pow p _ 3
    have hc : valAdd p (Units.mk0 (x ^ 3) hx3ne) = valAdd p (Units.mk0 x hx ^ 3) :=
      valAdd_eq_of_valuation_eq p _ _ (by simp)
    rw [hc, hmul]; norm_num
  -- 残りの 3 項はまとめて `ValAtLeast p (2 a)`
  have hrest : ValAtLeast p (2 * a) (W.a₂ * x ^ 2 + W.a₄ * x + W.a₆) := by
    have hx2 : ValAtLeast p (2 * a) (x ^ 2) := by
      have := valAtLeast_mul hxA hxA
      have hsq : x * x = x ^ 2 := by ring
      rw [hsq] at this
      refine valAtLeast_mono ?_ this; omega
    refine valAtLeast_add (valAtLeast_add ?_ ?_) ?_
    · refine valAtLeast_mono ?_ (valAtLeast_mul h2 hx2); omega
    · refine valAtLeast_mono ?_ (valAtLeast_mul h4 hxA); omega
    · refine valAtLeast_mono ?_ h6; omega
  have hsplit : x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆
      = x ^ 3 + (W.a₂ * x ^ 2 + W.a₄ * x + W.a₆) := by ring
  have hlt : valAdd p (Units.mk0 (x ^ 3) hx3ne) < 2 * a := by omega
  have hne0 : x ^ 3 + (W.a₂ * x ^ 2 + W.a₄ * x + W.a₆) ≠ 0 :=
    add_ne_zero_of_valAdd_lt hx3ne hrest hlt
  have hne : x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆ ≠ 0 := by rw [hsplit]; exact hne0
  refine ⟨hne, ?_⟩
  rw [valAdd_congr p hne hne0 hsplit, valAdd_add_eq_of_lt hx3ne hne0 hrest hlt, hx3]

/-- ★★★★★★★★**`v(x) < 0` なら `y ≠ 0`**。 -/
theorem y_ne_zero_of_valAdd_x_neg {x y : L} (hx : x ≠ 0)
    (heq : W.toAffine.Equation x y) (hneg : valAdd p (Units.mk0 x hx) < 0) : y ≠ 0 := by
  obtain ⟨hne, _⟩ := valAdd_rhs_eq p W hx hneg
  intro hy
  apply hne
  rw [WeierstrassCurve.Affine.equation_iff] at heq
  rw [← heq, hy]
  ring

/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 深い点では `2 v(y) = 3 v(x)`**（第 1070）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆Weierstrass 方程式の両辺で主要項が `y²` と `x³` しかないことによる。
★これが還元の核の濾過（形式群）の第 1 段である。 -/
theorem two_valAdd_y_eq_three_valAdd_x {x y : L} (hx : x ≠ 0) (hy : y ≠ 0)
    (heq : W.toAffine.Equation x y) (hneg : valAdd p (Units.mk0 x hx) < 0) :
    2 * valAdd p (Units.mk0 y hy) = 3 * valAdd p (Units.mk0 x hx) := by
  obtain ⟨h1, h2, h3, h4, h6⟩ := valAtLeast_coeffs p W
  obtain ⟨hRne, hR⟩ := valAdd_rhs_eq p W hx hneg
  set a : ℤ := valAdd p (Units.mk0 x hx) with ha_def
  set b : ℤ := valAdd p (Units.mk0 y hy) with hb_def
  have hxA : ValAtLeast p a x := valAtLeast_unit p (Units.mk0 x hx)
  have hyA : ValAtLeast p b y := valAtLeast_unit p (Units.mk0 y hy)
  have hy2 : ValAtLeast p (2 * b) (y ^ 2) := by
    have h := valAtLeast_mul hyA hyA
    have hsq : y * y = y ^ 2 := by ring
    rw [hsq] at h
    refine valAtLeast_mono ?_ h; omega
  have hx2 : ValAtLeast p (2 * a) (x ^ 2) := by
    have h := valAtLeast_mul hxA hxA
    have hsq : x * x = x ^ 2 := by ring
    rw [hsq] at h
    refine valAtLeast_mono ?_ h; omega
  have hx3 : ValAtLeast p (3 * a) (x ^ 3) := by
    have h := valAtLeast_mul hxA (valAtLeast_mul hxA hxA)
    have hcu : x * (x * x) = x ^ 3 := by ring
    rw [hcu] at h
    refine valAtLeast_mono ?_ h; omega
  have hxy : ValAtLeast p (a + b) (W.a₁ * x * y) := by
    have h := valAtLeast_mul (valAtLeast_mul h1 hxA) hyA
    refine valAtLeast_mono ?_ h; omega
  have ha3y : ValAtLeast p b (W.a₃ * y) := by
    have h := valAtLeast_mul h3 hyA
    refine valAtLeast_mono ?_ h; omega
  rw [WeierstrassCurve.Affine.equation_iff] at heq
  -- (1) 左辺の `valAdd` は `3a` なので `min (2b) (a+b) b ≤ 3a`
  have hLne : y ^ 2 + W.a₁ * x * y + W.a₃ * y ≠ 0 := by rw [heq]; exact hRne
  have hL : ValAtLeast p (min (2 * b) (min (a + b) b)) (y ^ 2 + W.a₁ * x * y + W.a₃ * y) := by
    refine valAtLeast_add (valAtLeast_add ?_ ?_) ?_
    · refine valAtLeast_mono ?_ hy2; omega
    · refine valAtLeast_mono ?_ hxy; omega
    · refine valAtLeast_mono ?_ ha3y; omega
  have hL' : min (2 * b) (min (a + b) b) ≤ 3 * a := by
    have hstep := hL hLne
    rwa [valAdd_congr p hLne hRne heq, hR] at hstep
  -- (2) `y²` の側から `min (3a) (min (a+b) b) ≤ 2b`
  have hy2eq : y ^ 2 = x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆
      + -(W.a₁ * x * y) + -(W.a₃ * y) := by linear_combination heq
  have hRfull : ValAtLeast p (min (3 * a) (min (a + b) b))
      (x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆ + -(W.a₁ * x * y) + -(W.a₃ * y)) := by
    refine valAtLeast_add (valAtLeast_add (valAtLeast_add (valAtLeast_add
      (valAtLeast_add ?_ ?_) ?_) ?_) ?_) ?_
    · refine valAtLeast_mono ?_ hx3; omega
    · refine valAtLeast_mono ?_ (valAtLeast_mul h2 hx2); omega
    · refine valAtLeast_mono ?_ (valAtLeast_mul h4 hxA); omega
    · refine valAtLeast_mono ?_ h6; omega
    · refine valAtLeast_mono ?_ (valAtLeast_neg hxy); omega
    · refine valAtLeast_mono ?_ (valAtLeast_neg ha3y); omega
  have hy2ne : y ^ 2 ≠ 0 := pow_ne_zero _ hy
  have hy2val : valAdd p (Units.mk0 (y ^ 2) hy2ne) = 2 * b := by
    have hmul : valAdd p (Units.mk0 y hy ^ 2) = ((2 : ℕ) : ℤ) * b := valAdd_pow p _ 2
    have hc : valAdd p (Units.mk0 (y ^ 2) hy2ne) = valAdd p (Units.mk0 y hy ^ 2) :=
      valAdd_eq_of_valuation_eq p _ _ (by simp)
    rw [hc, hmul]; norm_num
  have hR' : min (3 * a) (min (a + b) b) ≤ 2 * b := by
    have hne2 : x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆
        + -(W.a₁ * x * y) + -(W.a₃ * y) ≠ 0 := by rw [← hy2eq]; exact hy2ne
    have hstep := hRfull hne2
    rwa [valAdd_congr p hne2 hy2ne hy2eq.symm, hy2val] at hstep
  omega

/-- ★★★★★★★★★★★★**深さ `m` の取り出し**（第 1070）。

☆`v(x) = −2m`・`v(y) = −3m`（`m ≥ 1`）。 -/
theorem exists_depth_of_valAdd_x_neg {x y : L} (hx : x ≠ 0) (hy : y ≠ 0)
    (heq : W.toAffine.Equation x y) (hneg : valAdd p (Units.mk0 x hx) < 0) :
    ∃ m : ℕ, 0 < m ∧ valAdd p (Units.mk0 x hx) = -2 * (m : ℤ) ∧
      valAdd p (Units.mk0 y hy) = -3 * (m : ℤ) := by
  have hkey := two_valAdd_y_eq_three_valAdd_x p W hx hy heq hneg
  set a : ℤ := valAdd p (Units.mk0 x hx) with ha_def
  set b : ℤ := valAdd p (Units.mk0 y hy) with hb_def
  have hdvd : (2 : ℤ) ∣ a := by omega
  obtain ⟨k, hk⟩ := hdvd
  refine ⟨(-k).toNat, ?_, ?_, ?_⟩ <;> omega

end Point

/-! ## ★出典の紐付け(`.src`) -/

def two_valAdd_y_eq_three_valAdd_x.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い点では 2 v(y) = 3 v(x)——還元の核の濾過の第 1 段。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_depth_of_valAdd_x_neg.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い点の深さ m——v(x) = −2m・v(y) = −3m。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAdd_rhs_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v(x) < 0 なら方程式の右辺の valAdd は 3 v(x)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def y_ne_zero_of_valAdd_x_neg.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v(x) < 0 なら y ≠ 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def ValAtLeast.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(零を許す valAdd ≥ n——超距離の帳簿を omega に載せる)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
