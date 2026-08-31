/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.FaltingsWitness
import ABC3.Found.GaloisRep.SemistableFin
import ABC3.Found.GaloisRep.Compositum
import ABC3.Found.GaloisRep.HtFaltBounds
import ABC3.Found.GaloisRep.HtJBound
import ABC3.Found.GaloisRep.ResChar
import ABC3.Found.GaloisRep.ThetaDiscr
import ABC3.Found.GaloisRep.NorthcottHtJ
import ABC3.Found.GenEll.Velu
import ABC3.Meta.Claim

/-!
# `EllModuliData` の対象 —— 半安定な楕円曲線の族（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17–p.23。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★これは何か

`ResearchPaper/ellmoduli-witness-status.json` の `designChoice` で決めた設計を
**Lean の型として固定する**。★これで次のセッションが設計を検討し直さずに済み、
型検査で設計の齟齬も早く見つかる。

## ★★★設計

    Curve    := SSCurve（数体 `L ⊆ ℂ` と、その上の半安定な楕円曲線）
    EllClass := ℂ（`j` 不変量）
    cls E    := E.j

★`EllClass` は**本物の商**でなければならない——`EllClass := Curve` にすると
`northcott`（高さと次数で抑えた類の有限性）が偽になる（同型な Weierstrass モデルが
無限にあるから）。

★★`Curve` を**半安定なものだけ**に制限するのが要点である。そうすると
`EllModuliData` の `torsionExt` 群（`torsionExt` / `cls_torsionExt` /
`degOfDefinition_torsionExt` / `semiStable_torsionExt` / `hasMultRed_torsionExt` /
`primeToLocalHeights_torsionExt`）が `torsionExt := id` で**すべて自明に埋まる**
——これが唯一の mathlib 欠落（Néron–Ogg–Shafarevich、半安定還元の判定）だった。

☆ただしその結果は「半安定な曲線に対する `Lemma 3.5` 等」であって原文の主張そのもの
ではない（原文は `L(E[3],E[5])` への基底変換で一般の場合を半安定へ帰着させる）。
★したがって `.src` は条つきである。

## ★本ファイルで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `SSCurve` | ★数体とその上の半安定な楕円曲線 |
| `SSCurve.j` | ★`j` 不変量を `ℂ` へ |
| `SSCurve.deg` / `SSCurve.htFalt` / `SSCurve.degInf` | ★界面の欄に渡す量 |
| `SSCurve.deg_pos` | ★`degOfDefinition_pos` の中身 |
| `SSCurve.degInf_nonneg` | ★`deg∞ ≥ 0` |
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField

/-! ## ★★★★★★半安定な楕円曲線 -/

/-- ★★★★★★**数体 `L` とその上の半安定な楕円曲線**。

★`EllModuliData` の `Curve` 欄に渡す対象である。 -/
structure SSCurve where
  /-- ★定義体——**`ℂ` の中間体として取る**。

  ☆`Found/GaloisRep/NorthcottHtJ.lean` の `finite_j_of_htFalt_le`（`§9-1005`、★無条件）や
  `Found/GaloisRep/Lemma37C.lean` の `finite_j_of_condA` は、族の定義体を
  `fld : P → IntermediateField ℚ ℂ` の形で受け取る。★★ここを抽象的な `Type` ＋ 埋め込みに
  すると**それらが使えなくなる**（2026-08-31 の在庫の測定、第 753）。 -/
  K : IntermediateField ℚ ℂ
  [isNF : NumberField K]
  /-- 曲線。 -/
  W : WeierstrassCurve K
  [isEll : W.IsElliptic]
  /-- ★全ての有限素点で半安定。 -/
  ss : ∀ p : HeightOneSpectrum (𝓞 K), SemistableAt p W

attribute [instance] SSCurve.isNF SSCurve.isEll

namespace SSCurve

/-- 定義体（型として）。 -/
abbrev fld (E : SSCurve) : Type := E.K

/-- `j` を `ℂ` へ送る埋め込み（中間体の包含）。 -/
noncomputable def emb (E : SSCurve) : E.fld →+* ℂ := E.K.val.toRingHom

@[simp] theorem emb_apply (E : SSCurve) (x : E.fld) : E.emb x = (x : ℂ) := rfl

/-- ★★★★`j` 不変量（`ℂ` の元として）——`cls` 欄。 -/
noncomputable def j (E : SSCurve) : ℂ := E.emb (E.W.j)

/-- ★★★定義体の次数——`degOfDefinition` 欄。 -/
noncomputable def deg (E : SSCurve) : ℕ := Module.finrank ℚ E.fld

/-- ★★★★Faltings 高さ——`faltingsHeight` を作る材料。 -/
noncomputable def htFalt (E : SSCurve) : ℝ := htFaltOf E.fld E.W

/-- ★★★`deg_∞`——`degInf` を作る材料。 -/
noncomputable def degInf (E : SSCurve) : ℝ := degInfOf E.fld E.W

/-! ## ★★★★基本の性質 -/

/-- ★★`degOfDefinition_pos` 欄の中身。 -/
theorem deg_pos (E : SSCurve) : 0 < E.deg := Module.finrank_pos

/-- ★★`degInf_nonneg` 欄の中身。 -/
theorem degInf_nonneg (E : SSCurve) : 0 ≤ E.degInf := degInfOf_nonneg E.W

/-- ★★★★**変数変換で `ht^Falt` は変わらない**（第 722 を `SSCurve` の言葉で）。 -/
theorem htFalt_variableChange (E : SSCurve) (C : VariableChange E.fld)
    (hss : ∀ p : HeightOneSpectrum (𝓞 E.fld), SemistableAt p (C • E.W))
    (hell : (C • E.W).IsElliptic) :
    htFalt { K := E.K, W := C • E.W, isEll := hell, ss := hss }
      = E.htFalt := by
  show htFaltOf E.fld (C • E.W) = htFaltOf E.fld E.W
  exact htFaltOf_variableChange E.W C

end SSCurve

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★`ht^Falt` は `j` だけで決まる -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`j` が同じ `SSCurve` は `ht^Falt` も同じ**——★**無条件**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★これが `EllModuliData` の `faltingsHeight : EllClass → ℝ` の well-defined 性である
（`§9-1168`、第 741 の `htFaltOf_congr_j_of_emb` をそのまま `SSCurve` の言葉で）。 -/
theorem htFalt_congr_j (E E' : SSCurve) (hj : E.j = E'.j) : E.htFalt = E'.htFalt :=
  ABC3.Found.GaloisRep.htFaltOf_congr_j_of_emb E.fld E'.fld E.emb E'.emb E.W E'.W E.ss E'.ss hj

open scoped Classical in
/-- ★★★★★★★★★★★★★★**`EllModuliData` の `faltingsHeight` 欄**——`j` の函数として。

★`j` を持つ半安定曲線が無い `j` では `0` と定める（界面はそこでは何も要求しない）。 -/
noncomputable def faltingsHeightJ (j : ℂ) : ℝ :=
  if h : ∃ E : SSCurve, E.j = j then h.choose.htFalt else 0

/-- ★★★★★★★★★★★★★★★★★★**欄の値は曲線の `ht^Falt` に一致する**。 -/
theorem faltingsHeightJ_eq (E : SSCurve) : faltingsHeightJ E.j = E.htFalt := by
  classical
  have h : ∃ E' : SSCurve, E'.j = E.j := ⟨E, rfl⟩
  rw [faltingsHeightJ, dif_pos h]
  exact htFalt_congr_j h.choose E h.choose_spec

/-! ## ★★★★★★★★★★★★★★★★★★`deg∞` も `j` だけで決まる -/

/-- ★★★★★★★★★★★★★★**`j` が同じ `SSCurve` は `deg∞` も同じ**——★**無条件**。 -/
theorem degInf_congr_j (E E' : SSCurve) (hj : E.j = E'.j) : E.degInf = E'.degInf :=
  ABC3.Found.GaloisRep.degInfOf_congr_j_of_emb E.fld E'.fld E.emb E'.emb E.W E'.W E.ss E'.ss hj

open scoped Classical in
/-- ★★★★★★★★★★**`EllModuliData` の `degInf` 欄**——`j` の函数として。 -/
noncomputable def degInfJ (j : ℂ) : ℝ :=
  if h : ∃ E : SSCurve, E.j = j then h.choose.degInf else 0

/-- ★★★★★★★★★★★★**欄の値は曲線の `deg∞` に一致する**。 -/
theorem degInfJ_eq (E : SSCurve) : degInfJ E.j = E.degInf := by
  classical
  have h : ∃ E' : SSCurve, E'.j = E.j := ⟨E, rfl⟩
  rw [degInfJ, dif_pos h]
  exact degInf_congr_j h.choose E h.choose_spec

/-- ★★★★`deg∞ ≥ 0`——`degInf_nonneg` 欄。 -/
theorem degInfJ_nonneg (j : ℂ) : 0 ≤ degInfJ j := by
  classical
  by_cases h : ∃ E : SSCurve, E.j = j
  · rw [degInfJ, dif_pos h]
    exact h.choose.degInf_nonneg
  · rw [degInfJ, dif_neg h]

/-! ## ★★★★★★★★★★★★★★★★界面の 3 つの評価（類の水準で） -/

/-- ★★★★★★★★★★**`faltingsHeight` は下に有界**——`faltingsHeight_bddBelow` 欄。 -/
theorem faltingsHeightJ_bddBelow : ∃ B : ℝ, ∀ j : ℂ, B ≤ faltingsHeightJ j := by
  classical
  obtain ⟨B, hB⟩ := ABC3.Found.GaloisRep.exists_htFalt_bddBelow
  refine ⟨min B 0, fun j => ?_⟩
  by_cases h : ∃ E : SSCurve, E.j = j
  · rw [faltingsHeightJ, dif_pos h]
    exact le_trans (min_le_left _ _) (hB h.choose.fld h.choose.W)
  · rw [faltingsHeightJ, dif_neg h]
    exact min_le_right _ _

/-- ★★★★★★★★★★★★★★★★**`deg∞ ≤ 12·ht^Falt + A`**（類の水準で）——★**無条件**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★`htInf := 12·faltingsHeight` と取れば、これがそのまま `degInf_le_htInf` 欄である
（そして `htInf_bdeq_faltings` 欄は `C = 0` で自明になる）。 -/
theorem degInfJ_le_faltingsHeightJ :
    ∃ A : ℝ, 0 ≤ A ∧ ∀ j : ℂ, degInfJ j ≤ 12 * faltingsHeightJ j + A := by
  classical
  obtain ⟨A, hA0, hA⟩ := ABC3.Found.GaloisRep.exists_degInfOf_le_htFalt
  refine ⟨A, hA0, fun j => ?_⟩
  by_cases h : ∃ E : SSCurve, E.j = j
  · rw [degInfJ, dif_pos h, faltingsHeightJ, dif_pos h]
    exact hA h.choose.fld h.choose.W
  · rw [degInfJ, dif_neg h, faltingsHeightJ, dif_neg h]
    linarith

/-- ★★★★★★★★**`degInf_le_htInf` 欄の形**（`htInf := 12·faltingsHeight`）。 -/
theorem degInfJ_sub_htInfJ_le :
    ∃ C : ℝ, ∀ j : ℂ, degInfJ j - 12 * faltingsHeightJ j ≤ C := by
  obtain ⟨A, _, hA⟩ := degInfJ_le_faltingsHeightJ
  exact ⟨A, fun j => by linarith [hA j]⟩

/-! ## ★★★★★★★★★★★★`degLe` 欄と Northcott -/

/-- ★★★★**`EllModuliData` の `degLe` 欄**——`M_ell(ℚ̄)^{≤d}`。 -/
def degLeJ (d : ℕ) : Set ℂ := {x : ℂ | ∃ E : SSCurve, E.j = x ∧ E.deg ≤ d}

theorem mem_degLeJ (E : SSCurve) : E.j ∈ degLeJ E.deg := ⟨E, rfl, le_rfl⟩

/-! ## ★★★★★★★★★★★★述語の欄と `torsionExt` 群 -/

namespace SSCurve

/-- ★★★**素点 `p` での局所高さ**——半安定なので `v_p(Δ_min) = v_p(q)` である。 -/
noncomputable def localHeightAt (E : SSCurve) (p : HeightOneSpectrum (𝓞 E.fld)) : ℤ :=
  minDeltaExp p E.W

theorem localHeightAt_nonneg (E : SSCurve) (p : HeightOneSpectrum (𝓞 E.fld)) :
    0 ≤ E.localHeightAt p := minDeltaExp_nonneg p E.W

/-- ★★★★**`SemiStable` 欄**——`SSCurve` は構成上つねに半安定である。 -/
def SemiStable (E : SSCurve) : Prop := ∀ p : HeightOneSpectrum (𝓞 E.fld), SemistableAt p E.W

/-- ★★`SemiStable` 欄は無条件に成り立つ。 -/
theorem semiStable_all (E : SSCurve) : E.SemiStable := E.ss

/-- ★★★★**`HasMultRed` 欄**——半安定な曲線では「悪い還元」＝「乗法還元」である。 -/
def HasMultRed (E : SSCurve) : Prop := ∃ p : HeightOneSpectrum (𝓞 E.fld), E.localHeightAt p ≠ 0

/-- ★★★**`HasPotMultRed` 欄**——半安定な曲線では潜在的乗法還元は乗法還元と一致する
（変数変換で還元型は変わらず、半安定なら加法還元は起きない）。 -/
def HasPotMultRed (E : SSCurve) : Prop := E.HasMultRed

/-- ★★★★**`PrimeToLocalHeights` 欄**——`l` はどの局所高さとも素。 -/
def PrimeToLocalHeights (E : SSCurve) (l : ℕ) : Prop :=
  ∀ p : HeightOneSpectrum (𝓞 E.fld), E.localHeightAt p ≠ 0 →
    Nat.Coprime l (E.localHeightAt p).toNat

/-- ★★★★★★★★★★★★★★**`primeToLocalHeights_of_lt` 欄**——★**無条件**。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★中身は `HtFaltBounds.lean` の `coprime_minDeltaExp`（`§9-1052`）。 -/
theorem primeToLocalHeights_of_lt (E : SSCurve) (l : ℕ) (hl : Nat.Prime l)
    (h : (E.deg : ℝ) * degInfJ E.j < (l : ℝ) * Real.log 2) :
    E.PrimeToLocalHeights l := by
  rw [degInfJ_eq] at h
  exact fun p hp => coprime_minDeltaExp E.W E.W.isUnit_Δ.ne_zero l hl h p hp

end SSCurve

/-! ## ★★★★★★★★★★★★★★★★`torsionExt` 群——半安定に制限したので自明に埋まる -/

/-- ★★★★★★**`torsionExt` 欄**——`SSCurve` はすでに半安定なので**恒等写像**でよい。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★原文 p.20 の `L′`（3･5 捧れを有理化する次数 `23040` の拡大）は
「`E_{L′}` が全ての有限素点で半安定になるように」取るためのものである。
★★`Curve := SSCurve`（半安定なものだけ）と決めた時点で、この段は不要になる
——それが `ResearchPaper/ellmoduli-witness-status.json` の `semistableRestriction` の内容である。
☆代償として結論は「半安定な曲線について」に限定される。 -/
def torsionExt (E : SSCurve) : SSCurve := E

@[simp] theorem torsionExt_eq (E : SSCurve) : torsionExt E = E := rfl

theorem cls_torsionExt (E : SSCurve) : (torsionExt E).j = E.j := rfl

theorem degOfDefinition_torsionExt (E : SSCurve) :
    (torsionExt E).deg ≤ 23040 * E.deg := by
  show E.deg ≤ 23040 * E.deg
  have h := E.deg_pos
  omega

theorem semiStable_torsionExt (E : SSCurve) : (torsionExt E).SemiStable := E.semiStable_all

theorem hasMultRed_torsionExt (E : SSCurve) (h : E.HasPotMultRed) :
    (torsionExt E).HasMultRed := h

theorem primeToLocalHeights_torsionExt (E : SSCurve) (l : ℕ)
    (h : E.PrimeToLocalHeights l) (_ : Nat.Coprime l 30) :
    (torsionExt E).PrimeToLocalHeights l := h

/-! ## ★★★★★★★★★★★★★★★★`Curve` 欄は乗法還元を持つものに限る -/

/-- ★★★★★★★★★★**乗法還元を持つ半安定曲線**。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★★★`Check/GenEll/EllModuliDegInfPos.lean` の測定（2026-08-31、第 745）:
`EllModuliData` の `multCard_pos`・`localHt_pos`・`multPrime_prime`・`sum_localHt_eq` は
**どの `E : Curve` に対しても `deg∞(cls E) > 0` を強制する**。
★したがって `Curve` 欄は「至る所良還元」の曲線を含めない——原文の
*Degenerating* Elliptic Curves という題そのものである。 -/
structure DegCurve where
  /-- 半安定な楕円曲線。 -/
  toSSCurve : SSCurve
  /-- ★少なくとも 1 つの素点で乗法還元。 -/
  multRed : toSSCurve.HasMultRed

namespace DegCurve

/-- `j` 不変量。 -/
noncomputable def j (E : DegCurve) : ℂ := E.toSSCurve.j

/-- 定義体の次数。 -/
noncomputable def deg (E : DegCurve) : ℕ := E.toSSCurve.deg

/-- ★★★★★★**`deg∞ > 0`**——界面が強制する条件を満たす。 -/
theorem degInf_pos (E : DegCurve) : 0 < degInfJ E.j := by
  obtain ⟨p, hp⟩ := E.multRed
  rw [DegCurve.j, degInfJ_eq]
  exact degInfOf_pos_of_minDeltaExp_ne_zero E.toSSCurve.W
    E.toSSCurve.W.isUnit_Δ.ne_zero p hp

theorem deg_pos (E : DegCurve) : 0 < E.deg := E.toSSCurve.deg_pos

end DegCurve

/-! ## ★★★★★★★★★★★★★★★★§4 の帳簿 -/

namespace SSCurve

open scoped Classical in
/-- ★★★**悪い還元の素点の有限集合**。 -/
noncomputable def badFinset (E : SSCurve) : Finset (HeightOneSpectrum (𝓞 E.fld)) :=
  (minDeltaExp_finite E.W E.W.isUnit_Δ.ne_zero).toFinset

open scoped Classical in
theorem mem_badFinset (E : SSCurve) (p : HeightOneSpectrum (𝓞 E.fld)) :
    p ∈ E.badFinset ↔ minDeltaExp p E.W ≠ 0 := by
  rw [badFinset, Set.Finite.mem_toFinset]
  rfl

open scoped Classical in
/-- ★★★★★★★**帳簿の主等式**——`∑_{p bad} v_p(Δ_min)·log N(p) = d·deg∞`。 -/
theorem sum_badFinset_eq (E : SSCurve) :
    ∑ p ∈ E.badFinset, (minDeltaExp p E.W : ℝ) * Real.log (Ideal.absNorm p.asIdeal)
      = (E.deg : ℝ) * degInfJ E.j := by
  rw [degInfJ_eq]
  show _ = (Module.finrank ℚ E.fld : ℝ) * degInfOf E.fld E.W
  rw [finrank_mul_degInfOf]
  refine (finsum_eq_finset_sum_of_support_subset _ ?_).symm
  intro q hq
  simp only [Function.mem_support, ne_eq, mul_eq_zero, not_or] at hq
  refine Finset.mem_coe.2 ((mem_badFinset E q).2 ?_)
  exact_mod_cast hq.1

end SSCurve

namespace DegCurve

open scoped Classical in
/-- ★★★**`multCard` 欄**——悪い還元の素点の個数。 -/
noncomputable def multCard (E : DegCurve) : ℕ := E.toSSCurve.badFinset.card

open scoped Classical in
/-- ★★★★**`multCard_pos` 欄**——乗法還元を持つので空でない。 -/
theorem multCard_pos (E : DegCurve) : 0 < E.multCard := by
  obtain ⟨p, hp⟩ := E.multRed
  exact Finset.card_pos.2 ⟨p, (SSCurve.mem_badFinset _ p).2 hp⟩

open scoped Classical in
/-- ★★悪い素点の番号付け。 -/
noncomputable def multIdx (E : DegCurve) (j : Fin E.multCard) :
    HeightOneSpectrum (𝓞 E.toSSCurve.fld) :=
  (E.toSSCurve.badFinset.equivFin.symm j : _)

open scoped Classical in
theorem multIdx_mem (E : DegCurve) (j : Fin E.multCard) :
    E.multIdx j ∈ E.toSSCurve.badFinset := (E.toSSCurve.badFinset.equivFin.symm j).2

open scoped Classical in
/-- ★★★**`multPrime` 欄**——その素点の剰余標数。 -/
noncomputable def multPrime (E : DegCurve) (j : Fin E.multCard) : ℕ :=
  resChar (E.multIdx j)

open scoped Classical in
theorem multPrime_prime (E : DegCurve) (j : Fin E.multCard) : Nat.Prime (E.multPrime j) :=
  resChar_prime _

open scoped Classical in
/-- ★★★★**`localHt` 欄**——`23040·v_p(Δ_min)·f_p`。

★係数 `23040` は原文 p.23 の『the “`h`” of Lemma 4.2 corresponds to
`23040d · deg∞([E_L])`』から来る。 -/
noncomputable def localHt (E : DegCurve) (j : Fin E.multCard) : ℕ :=
  23040 * (minDeltaExp (E.multIdx j) E.toSSCurve.W).toNat
    * ((Ideal.span {(resChar (E.multIdx j) : ℤ)}).inertiaDeg (E.multIdx j).asIdeal)

open scoped Classical in
theorem localHt_pos (E : DegCurve) (j : Fin E.multCard) : 0 < E.localHt j := by
  have h1 : minDeltaExp (E.multIdx j) E.toSSCurve.W ≠ 0 :=
    (SSCurve.mem_badFinset _ _).1 (E.multIdx_mem j)
  have h2 : 0 < (minDeltaExp (E.multIdx j) E.toSSCurve.W).toNat := by
    have := minDeltaExp_nonneg (E.multIdx j) E.toSSCurve.W
    omega
  have h3 := inertiaDeg_pos (E.multIdx j)
  rw [localHt]
  positivity

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★
**`sum_localHt_eq` 欄**——`∑ h_j·log(p_j) = 23040·d·deg∞([E_L])`（原文 p.23）。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★機構は `log N(p) = f_p·log(剰余標数)`（`ResChar.lean`、`§9-1173`）と
`∑_{p bad} v_p·log N(p) = d·deg∞`（`sum_badFinset_eq`）。 -/
theorem sum_localHt_eq (E : DegCurve) :
    (∑ j : Fin E.multCard, (E.localHt j : ℝ) * Real.log (E.multPrime j))
      = 23040 * (E.deg : ℝ) * degInfJ E.j := by
  have hterm : ∀ j : Fin E.multCard, (E.localHt j : ℝ) * Real.log (E.multPrime j)
      = 23040 * ((minDeltaExp (E.multIdx j) E.toSSCurve.W : ℝ)
          * Real.log (Ideal.absNorm (E.multIdx j).asIdeal)) := by
    intro j
    have hnn := minDeltaExp_nonneg (E.multIdx j) E.toSSCurve.W
    have hcast : (((minDeltaExp (E.multIdx j) E.toSSCurve.W).toNat : ℕ) : ℝ)
        = (minDeltaExp (E.multIdx j) E.toSSCurve.W : ℝ) := by
      exact_mod_cast congrArg (fun z : ℤ => (z : ℝ)) (Int.toNat_of_nonneg hnn)
    rw [localHt, multPrime, log_absNorm_eq (E.multIdx j)]
    push_cast
    rw [hcast]
    ring
  rw [Finset.sum_congr rfl (fun j _ => hterm j), ← Finset.mul_sum]
  have hsum : (∑ j : Fin E.multCard, (minDeltaExp (E.multIdx j) E.toSSCurve.W : ℝ)
      * Real.log (Ideal.absNorm (E.multIdx j).asIdeal))
      = ∑ p ∈ E.toSSCurve.badFinset,
          (minDeltaExp p E.toSSCurve.W : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
    rw [← Finset.sum_coe_sort E.toSSCurve.badFinset
      (fun p => (minDeltaExp p E.toSSCurve.W : ℝ) * Real.log (Ideal.absNorm p.asIdeal))]
    exact Equiv.sum_comp E.toSSCurve.badFinset.equivFin.symm
      (fun p : {x // x ∈ E.toSSCurve.badFinset} =>
        (minDeltaExp (p : HeightOneSpectrum (𝓞 E.toSSCurve.fld)) E.toSSCurve.W : ℝ)
          * Real.log (Ideal.absNorm (p : HeightOneSpectrum (𝓞 E.toSSCurve.fld)).asIdeal))
  rw [hsum, SSCurve.sum_badFinset_eq]
  show 23040 * ((E.toSSCurve.deg : ℝ) * degInfJ E.toSSCurve.j) = _
  rw [DegCurve.deg, DegCurve.j]
  ring

end DegCurve

/-! ## ★★★★★★★★★★★★★★★★★★`badPrimes` 欄 -/

section SumHelpers

variable {ι α : Type*} [DecidableEq α]

/-- ★非負なら合併の和は和の和以下。 -/
theorem sum_union_le_of_nonneg (s t : Finset α) (g : α → ℝ) (hg : ∀ a, 0 ≤ g a) :
    ∑ a ∈ s ∪ t, g a ≤ ∑ a ∈ s, g a + ∑ a ∈ t, g a := by
  have h := Finset.sum_union_inter (s₁ := s) (s₂ := t) (f := g)
  have h2 : 0 ≤ ∑ a ∈ s ∩ t, g a := Finset.sum_nonneg (fun a _ => hg a)
  linarith

/-- ★非負なら `biUnion` の和は和の和以下。 -/
theorem sum_biUnion_le_of_nonneg (s : Finset ι) (t : ι → Finset α) (g : α → ℝ)
    (hg : ∀ a, 0 ≤ g a) [DecidableEq ι] :
    ∑ a ∈ s.biUnion t, g a ≤ ∑ i ∈ s, ∑ a ∈ t i, g a := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert i s hi ih =>
    rw [Finset.biUnion_insert, Finset.sum_insert hi]
    have h1 := sum_union_le_of_nonneg (t i) (s.biUnion t) g hg
    linarith

/-- ★`insert` の和は 1 項＋残りの和以下（非負のとき）。 -/
theorem sum_insert_le_of_nonneg (a : α) (s : Finset α) (g : α → ℝ) (hg : ∀ x, 0 ≤ g x) :
    ∑ x ∈ insert a s, g x ≤ g a + ∑ x ∈ s, g x := by
  classical
  by_cases h : a ∈ s
  · rw [Finset.insert_eq_self.2 h]
    have : 0 ≤ g a := hg a
    linarith
  · rw [Finset.sum_insert h]

end SumHelpers

namespace DegCurve

open scoped Classical in
/-- ★★★★**`badPrimes` 欄**——乗法還元の素点の下の有理素数と、
局所高さの素因数分解に現れる素数。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★原文 p.22:『write `S∘` for the union of `S`, the primes of `ℚ` that lie under primes of
potentially multiplicative reduction of `E_L`, and the primes that appear in the
prime decomposition of the local heights of `E_L`』。 -/
noncomputable def badPrimes (E : DegCurve) : Finset ℕ :=
  Finset.univ.biUnion
    (fun j : Fin E.multCard => insert (E.multPrime j) (E.localHt j).primeFactors)

open scoped Classical in
theorem badPrimes_prime (E : DegCurve) : ∀ q ∈ E.badPrimes, Nat.Prime q := by
  intro q hq
  obtain ⟨j, -, hj⟩ := Finset.mem_biUnion.1 hq
  rcases Finset.mem_insert.1 hj with h | h
  · rw [h]; exact E.multPrime_prime j
  · exact Nat.prime_of_mem_primeFactors h

open scoped Classical in
/-- ★★★★**`PrimeToMultPrimes` 欄**。 -/
def PrimeToMultPrimes (E : DegCurve) (l : ℕ) : Prop := ∀ j : Fin E.multCard, l ≠ E.multPrime j

open scoped Classical in
/-- ★★悪い素点はすべて番号付けに現れる。 -/
theorem exists_multIdx (E : DegCurve) (p : HeightOneSpectrum (𝓞 E.toSSCurve.fld))
    (hp : p ∈ E.toSSCurve.badFinset) : ∃ j : Fin E.multCard, E.multIdx j = p := by
  refine ⟨E.toSSCurve.badFinset.equivFin ⟨p, hp⟩, ?_⟩
  show ((E.toSSCurve.badFinset.equivFin.symm
    (E.toSSCurve.badFinset.equivFin ⟨p, hp⟩) : _) : _) = p
  rw [Equiv.symm_apply_apply]

open scoped Classical in
/-- ★★★★★★★★★★★★**`primeTo_badPrimes` 欄**——★**無条件**。 -/
theorem primeTo_badPrimes (E : DegCurve) (l : ℕ) (hl : Nat.Prime l) (h : l ∉ E.badPrimes) :
    E.PrimeToMultPrimes l ∧ E.toSSCurve.PrimeToLocalHeights l := by
  constructor
  · intro j hj
    exact h (Finset.mem_biUnion.2 ⟨j, Finset.mem_univ j,
      Finset.mem_insert.2 (Or.inl hj)⟩)
  · intro p hp
    obtain ⟨j, hj⟩ := exists_multIdx E p ((SSCurve.mem_badFinset _ p).2 hp)
    have hdvd : (E.toSSCurve.localHeightAt p).toNat ∣ E.localHt j := by
      rw [localHt, hj]
      show (minDeltaExp p E.toSSCurve.W).toNat ∣ _
      exact ⟨23040 * ((Ideal.span {(resChar p : ℤ)}).inertiaDeg p.asIdeal), by ring⟩
    rw [Nat.Prime.coprime_iff_not_dvd hl]
    intro hcon
    exact h (Finset.mem_biUnion.2 ⟨j, Finset.mem_univ j, Finset.mem_insert.2 (Or.inr
      (Nat.mem_primeFactors.2 ⟨hl, dvd_trans hcon hdvd, (E.localHt_pos j).ne'⟩))⟩)

open scoped Classical in
/-- ★★★★★★★★★★★★★★**`sum_log_badPrimes_le` 欄**——★**無条件**。

★`∏_{p ∣ n} p ∣ n`（`Nat.prod_primeFactors_dvd`）から
`∑_{p ∣ n} log p ≤ log n ≤ log(n+1)`。 -/
theorem sum_log_badPrimes_le (E : DegCurve) :
    (∑ q ∈ E.badPrimes, Real.log q)
      ≤ (∑ j : Fin E.multCard, Real.log (E.multPrime j))
        + (∑ j : Fin E.multCard, Real.log ((E.localHt j : ℝ) + 1)) := by
  have hnn : ∀ q : ℕ, 0 ≤ Real.log q := by
    intro q
    rcases Nat.eq_zero_or_pos q with h | h
    · rw [h]; simp
    · exact Real.log_nonneg (by exact_mod_cast h)
  have hstep : ∀ j : Fin E.multCard,
      (∑ q ∈ insert (E.multPrime j) (E.localHt j).primeFactors, Real.log q)
        ≤ Real.log (E.multPrime j) + Real.log ((E.localHt j : ℝ) + 1) := by
    intro j
    refine le_trans (sum_insert_le_of_nonneg _ _ _ hnn) ?_
    have hpf : (∑ q ∈ (E.localHt j).primeFactors, Real.log q)
        ≤ Real.log ((E.localHt j : ℝ) + 1) := by
      have hprod : ∏ q ∈ (E.localHt j).primeFactors, q ∣ E.localHt j :=
        Nat.prod_primeFactors_dvd _
      have hle : (∏ q ∈ (E.localHt j).primeFactors, q) ≤ E.localHt j :=
        Nat.le_of_dvd (E.localHt_pos j) hprod
      have hlog : (∑ q ∈ (E.localHt j).primeFactors, Real.log q)
          = Real.log ((∏ q ∈ (E.localHt j).primeFactors, q : ℕ) : ℝ) := by
        push_cast
        rw [Real.log_prod]
        intro q hq
        exact_mod_cast (Nat.prime_of_mem_primeFactors hq).ne_zero
      rw [hlog]
      refine Real.log_le_log ?_ ?_
      · have : 0 < ∏ q ∈ (E.localHt j).primeFactors, q :=
          Finset.prod_pos (fun q hq => (Nat.prime_of_mem_primeFactors hq).pos)
        exact_mod_cast this
      · have : ((∏ q ∈ (E.localHt j).primeFactors, q : ℕ) : ℝ) ≤ (E.localHt j : ℝ) := by
          exact_mod_cast hle
        linarith
    linarith
  refine le_trans (sum_biUnion_le_of_nonneg _ _ _ hnn) ?_
  rw [← Finset.sum_add_distrib]
  exact Finset.sum_le_sum (fun j _ => hstep j)

end DegCurve

/-! ## ★★★★★★★★★★★★`GaloisFinite` と `CompactlyBounded` -/

/-- ★★★★**`GaloisFinite` 欄**——類の集合が有限であること。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

## ★★★★★2026-08-31 の訂正（第 756）

☆以前ここは `GaloisFiniteJ S := S.Finite` としていたが、**それでは `lcyclicExc` が作れない**。
`lcyclicExc` は「`ht^Falt` が有界な類の集合」であり、**次数を止めなければ無限集合**である。

★原文の `Lemma 3.7` の例外集合は `Exc_d`——**次数 `d` を止めた形**である
（`Found/GaloisRep/Lemma37C.lean` の `finite_j_of_condA` も `d` を受け取る）。
★★そこで「各 `d` について `S ∩ M_ell(ℚ̄)^{≤d}` が有限」を `GaloisFinite` の内容とする。

☆原文の Galois-finite（有限個の Galois 軌道の合併）はこれより**強い**——逸脱として記録する。 -/
def GaloisFiniteJ (S : Set ℂ) : Prop := ∀ d : ℕ, (S ∩ degLeJ d).Finite

/-- ★★**`galoisFinite_union` 欄**。 -/
theorem galoisFiniteJ_union (S T : Set ℂ) (hS : GaloisFiniteJ S) (hT : GaloisFiniteJ T) :
    GaloisFiniteJ (S ∪ T) := by
  intro d
  rw [Set.union_inter_distrib_right]
  exact (hS d).union (hT d)

theorem galoisFiniteJ_empty : GaloisFiniteJ (∅ : Set ℂ) := by
  intro d
  simp

/-- ★★★★**`CompactlyBounded` 欄**——アルキメデス側の高さが有界であること。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★原文の compactly bounded（`Example 1.3, (ii)`）は有限個の素点 `V` と
そこでのコンパクト集合 `K_v` で押さえることである。
★★`EllClass := ℂ` の下でアルキメデス素点の条件が言っているのは
「`j` の**すべての共役**が有界」——すなわち `h_∞(j)` が有界——ということなので、
それを直接述語にする。

☆`j` そのものの有界性（`Bornology.IsBounded`）では**足りない**
——共役は有界にならないからである。 -/
def CompactlyBoundedJ (S : Set ℂ) : Prop :=
  ∃ M : ℝ, ∀ E : SSCurve, E.j ∈ S → htArchJ E.fld E.W ≤ M

/-- ★**`compactlyBounded_empty` 欄**。 -/
theorem compactlyBoundedJ_empty : CompactlyBoundedJ (∅ : Set ℂ) :=
  ⟨0, fun _ h => absurd h (Set.notMem_empty _)⟩

/-! ## ★★★★★★★★★★★★★★`noMultRedExc` 群——`DegCurve` なら自明 -/

/-- ★★★★**`noMultRedExc` 欄**——空集合でよい。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★★`Curve := DegCurve`（乗法還元を持つ半安定曲線）と決めたので、
`¬ HasMultRed E` を満たす `E : Curve` は**存在しない**。
★したがって `mem_noMultRedExc` 欄は空虚に成り立ち、例外集合は空でよい。
☆これは `Check/GenEll/EllModuliDegInfPos.lean` の測定（界面が `deg∞ > 0` を強制する）
から強いられた `Curve` の制限の**副産物**である。 -/
def noMultRedExcJ (_ : Set ℂ) : Set ℂ := ∅

theorem galoisFiniteJ_noMultRedExcJ (KV : Set ℂ) (_ : CompactlyBoundedJ KV) :
    GaloisFiniteJ (noMultRedExcJ KV) := galoisFiniteJ_empty

theorem mem_noMultRedExcJ (KV : Set ℂ) (E : DegCurve) (_ : E.j ∈ KV)
    (h : ¬ E.toSSCurve.HasMultRed) : E.j ∈ noMultRedExcJ KV :=
  absurd E.multRed h

/-! ## ★★★★★★`MinimalField` 欄 -/

/-- ★★★**`MinimalField` 欄**——定義体の次数が最小であること。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★原文『`L` is a minimal field of definition of the point `[E_L]`』。
☆`Skeleton/GenEll/Section4.lean` の測定（第 88 行）によれば
`Corollary 4.3` の証明はこの仮説を**使っていない**。 -/
def MinimalFieldJ (E : DegCurve) : Prop := ∀ E' : DegCurve, E'.j = E.j → E.deg ≤ E'.deg

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`northcott` 欄**——★**無条件**。

原文 (GenEll p.17):
> on Mell(Q). In particular, if C ∈R, then the set of points [E] ∈Mell(Q)≤d such

★★★在庫の測定（2026-08-31、第 753）: `Found/GaloisRep/NorthcottHtJ.lean` の
`finite_j_of_htFalt_le`（`§9-1005`）が**すでに無条件で証明されている**。
☆それに気づかず `Skeleton/GenEll/NorthcottJ.lean` に節点を立てていた（第 743）が、
**不要だった**——本定理がその節点を消す。

☆`finite_j_of_htFalt_le` は族の定義体を `fld : P → IntermediateField ℚ ℂ` の形で受ける。
★これが `SSCurve` の `K : IntermediateField ℚ ℂ` という設計の理由である。 -/
theorem northcottJ (C : ℝ) (d : ℕ) :
    {x ∈ degLeJ d | faltingsHeightJ x ≤ C}.Finite := by
  classical
  refine Set.Finite.subset
    (ABC3.Found.GaloisRep.finite_j_of_htFalt_le (P := {E : SSCurve // E.deg ≤ d}) d
      (fun E => E.1.K) (fun E => E.1.isNF) (fun E => E.2)
      (fun E => E.1.W) (fun E => E.1.isEll)
      (fun E => E.1.htFalt) (fun _ => rfl) C) ?_
  rintro x ⟨⟨E, hEj, hEd⟩, hC⟩
  refine ⟨⟨E, hEd⟩, ?_, ?_⟩
  · show E.htFalt ≤ C
    rw [← faltingsHeightJ_eq E, hEj]
    exact hC
  · show ((E.W.j : E.K) : ℂ) = x
    rw [← hEj]
    rfl


/-- ★★★★★★★★**`ht^Falt` が有界な集合は `GaloisFinite`**——★**無条件**。

★`northcottJ`（`§9-1179`、第 753）そのものである。
★★これが `galoisFinite_lcyclicExc` 欄の中身になる。 -/
theorem galoisFiniteJ_htFalt_le (B : ℝ) : GaloisFiniteJ {x : ℂ | faltingsHeightJ x ≤ B} := by
  intro d
  refine (northcottJ B d).subset ?_
  rintro x ⟨hB, hd⟩
  exact ⟨hd, hB⟩

/-! ## ★★★★★★★★★★★★★★★★★★★★★★`Curve` 欄は「実現される類」である -/

open scoped Classical in
/-- ★★★★★**`j` を実現する `DegCurve` のうち定義体の次数が最小のもの**。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★★原文 `Corollary 4.3` の仮説『`L` is a **minimal field of definition** of the point
`[E_L]`』そのものである。 -/
noncomputable def repCurve {x : ℂ} (h : ∃ E : DegCurve, E.j = x) : DegCurve :=
  (Nat.find_spec (p := fun n => ∃ E : DegCurve, E.j = x ∧ E.deg = n)
    ⟨h.choose.deg, h.choose, h.choose_spec, rfl⟩).choose

open scoped Classical in
theorem repCurve_j {x : ℂ} (h : ∃ E : DegCurve, E.j = x) : (repCurve h).j = x :=
  (Nat.find_spec (p := fun n => ∃ E : DegCurve, E.j = x ∧ E.deg = n)
    ⟨h.choose.deg, h.choose, h.choose_spec, rfl⟩).choose_spec.1

open scoped Classical in
theorem repCurve_deg {x : ℂ} (h : ∃ E : DegCurve, E.j = x) :
    (repCurve h).deg
      = Nat.find (p := fun n => ∃ E : DegCurve, E.j = x ∧ E.deg = n)
          ⟨h.choose.deg, h.choose, h.choose_spec, rfl⟩ :=
  (Nat.find_spec (p := fun n => ∃ E : DegCurve, E.j = x ∧ E.deg = n)
    ⟨h.choose.deg, h.choose, h.choose_spec, rfl⟩).choose_spec.2

open scoped Classical in
/-- ★★★★★**代表元の定義体の次数は最小である**。 -/
theorem repCurve_deg_le {x : ℂ} (h : ∃ E : DegCurve, E.j = x) (E : DegCurve) (hE : E.j = x) :
    (repCurve h).deg ≤ E.deg := by
  rw [repCurve_deg h]
  exact Nat.find_le ⟨E, hE, rfl⟩

/-- ★★★★★★★★★★★★★★★★**`EllModuliData` の `Curve` 欄**——実現される類。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

## ★★★★★なぜ「曲線」ではなく「類」を `Curve` に取るか

`EllModuliData` の `logDiffMell : EllClass → ℝ` は**類の函数**でなければならないが、
`log-diff` は**定義体の分岐**で決まる量である。★同じ `j` を持つ曲線でも定義体は
いくらでも大きく取れるので、`sum_log_ramPrimes_le` 欄（`∑_{ramPrimes} log p ≤
∑_{badPrimes} log p + 3d·log-diff`）は「勝手な定義体の曲線」に対しては**破れる**。

★★そこで `Curve` を**実現される類**とし、各類に**定義体の次数が最小の代表元**を
取ることにする。これで `degOfDefinition`・`logDiffMell` などがすべて類の函数になる。

★★★これは原文に忠実である——`Corollary 4.3` の仮説はまさに
『`L` is a **minimal field of definition** of the point `[E_L]`』である。 -/
def RealizedClass : Type := {x : ℂ // ∃ E : DegCurve, E.j = x}

namespace RealizedClass

/-- ★★★★代表元。 -/
noncomputable def rep (x : RealizedClass) : DegCurve := repCurve x.2

theorem rep_j (x : RealizedClass) : x.rep.j = x.1 := repCurve_j x.2

theorem rep_deg_le (x : RealizedClass) (E : DegCurve) (hE : E.j = x.1) :
    x.rep.deg ≤ E.deg := repCurve_deg_le x.2 E hE

/-- ★★**`cls` 欄**。 -/
def cls (x : RealizedClass) : ℂ := x.1

/-- ★★**`degOfDefinition` 欄**。 -/
noncomputable def degOfDefinition (x : RealizedClass) : ℕ := x.rep.deg

/-- ★**`degOfDefinition_pos` 欄**。 -/
theorem degOfDefinition_pos (x : RealizedClass) : 0 < x.degOfDefinition := x.rep.deg_pos

/-- ★★★**`MinimalField` 欄**——代表元は最小なので**つねに成り立つ**。 -/
theorem minimalField (x : RealizedClass) : MinimalFieldJ x.rep :=
  fun E' hE' => x.rep_deg_le E' (by rw [hE', x.rep_j])

/-- ★★★★**`deg∞` は代表元で測っても同じ**。 -/
theorem degInfJ_rep (x : RealizedClass) : degInfJ x.rep.j = degInfJ x.cls := by
  rw [x.rep_j]; rfl

/-- ★★★★★★★★★★**`sum_localHt_eq` 欄**（類の水準で）。 -/
theorem sum_localHt_eq (x : RealizedClass) :
    (∑ j : Fin x.rep.multCard, (x.rep.localHt j : ℝ) * Real.log (x.rep.multPrime j))
      = 23040 * (x.degOfDefinition : ℝ) * degInfJ x.cls := by
  rw [DegCurve.sum_localHt_eq x.rep, degOfDefinition, degInfJ_rep]

/-- ★★★★**`degInf` は正**（`Curve` の制限から）。 -/
theorem degInf_pos (x : RealizedClass) : 0 < degInfJ x.cls := by
  rw [← degInfJ_rep]
  exact x.rep.degInf_pos

end RealizedClass

/-! ## ★★★★★★★★★★★★★★★★★★★★`logDiffMell` と `ramPrimes` 群 -/

/-- ★★★素因数の対数の和は本体の対数以下（根基は本体を割るから）。 -/
theorem sum_log_primeFactors_le (n : ℕ) (hn : 0 < n) :
    (∑ q ∈ n.primeFactors, Real.log q) ≤ Real.log n := by
  have hprod : ∏ q ∈ n.primeFactors, q ∣ n := Nat.prod_primeFactors_dvd _
  have hle : (∏ q ∈ n.primeFactors, q) ≤ n := Nat.le_of_dvd hn hprod
  have hpos : 0 < ∏ q ∈ n.primeFactors, q :=
    Finset.prod_pos (fun q hq => (Nat.prime_of_mem_primeFactors hq).pos)
  have hlog : (∑ q ∈ n.primeFactors, Real.log q)
      = Real.log ((∏ q ∈ n.primeFactors, q : ℕ) : ℝ) := by
    push_cast
    rw [Real.log_prod]
    intro q hq
    exact_mod_cast (Nat.prime_of_mem_primeFactors hq).ne_zero
  rw [hlog]
  refine Real.log_le_log (by exact_mod_cast hpos) ?_
  exact_mod_cast hle

namespace RealizedClass

/-- ★★★★★★**`logDiffMell` 欄**——`log-diff = log|disc L| / d`。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★`Definition 1.5, (iii)` を `M̄_ell` に適用したもの。
★★`Curve` を**実現される類**に取った（`RealizedClass`）ので、
定義体は類ごとに一意に決まり、これは**類の函数**になる。 -/
noncomputable def logDiffMell (x : RealizedClass) : ℝ :=
  Real.log ((|NumberField.discr x.rep.toSSCurve.fld| : ℤ) : ℝ) / (x.degOfDefinition : ℝ)

theorem logDiffMell_nonneg (x : RealizedClass) : 0 ≤ x.logDiffMell := by
  refine div_nonneg ?_ (by positivity)
  refine Real.log_nonneg ?_
  have : (1:ℤ) ≤ |NumberField.discr x.rep.toSSCurve.fld| :=
    Int.one_le_abs (NumberField.discr_ne_zero (K := x.rep.toSSCurve.fld))
  exact_mod_cast this

open scoped Classical in
/-- ★★★★★**`ramPrimes` 欄**——`badPrimes` に、`L` で分岐する素数
（`disc L` を割る素数）と、分岐指数を割りうる素数（`d` 以下の素数）を加えたもの。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★原文 p.22:『write `S•` for the union of `S∘`, the primes of `ℚ` that **ramify in `L`**,
and the primes that divide the **ramification indices** of primes of `ℚ` in `L`』。
☆分岐指数は `≤ d` なので、それを割る素数は `d` 以下である。 -/
noncomputable def ramPrimes (x : RealizedClass) : Finset ℕ :=
  x.rep.badPrimes ∪ (NumberField.discr x.rep.toSSCurve.fld).natAbs.primeFactors
    ∪ ((Finset.range (x.degOfDefinition + 1)).filter Nat.Prime)

open scoped Classical in
theorem ramPrimes_prime (x : RealizedClass) : ∀ q ∈ x.ramPrimes, Nat.Prime q := by
  intro q hq
  rcases Finset.mem_union.1 hq with h | h
  · rcases Finset.mem_union.1 h with h' | h'
    · exact x.rep.badPrimes_prime q h'
    · exact Nat.prime_of_mem_primeFactors h'
  · exact (Finset.mem_filter.1 h).2

open scoped Classical in
theorem badPrimes_subset_ramPrimes (x : RealizedClass) : x.rep.badPrimes ⊆ x.ramPrimes :=
  fun _ hq => Finset.mem_union.2 (Or.inl (Finset.mem_union.2 (Or.inl hq)))

open scoped Classical in
/-- ★★★★**`PrimeToRamification` 欄**。 -/
def PrimeToRamification (x : RealizedClass) (l : ℕ) : Prop := l ∉ x.ramPrimes

open scoped Classical in
theorem primeTo_ramPrimes (x : RealizedClass) (l : ℕ) (_ : Nat.Prime l)
    (h : l ∉ x.ramPrimes) : x.PrimeToRamification l := h

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★
**`sum_log_ramPrimes_le` 欄**——★**無条件**。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★★係数 `3` の内訳:

* 分岐する素数の分 `∑_{p ∣ disc} log p ≤ log|disc|` —— **1 単位**
* `d` 以下の素数の分 `θ(d) ≤ 2·log|disc|`（`ThetaDiscr.lean`、`§9-1176`）—— **2 単位**

☆原文『since [as is easily verified, by considering the **trace** of an extension of
number fields] the primes appearing in the arithmetic divisor that gives rise to
“`log-diff_Mell`” appear with **multiplicity ≥ one less than the ramification indices**』
に対応する段を、Hermite–Minkowski と Chebyshev で置き換えた形である。 -/
theorem sum_log_ramPrimes_le (x : RealizedClass) :
    (∑ q ∈ x.ramPrimes, Real.log q)
      ≤ (∑ q ∈ x.rep.badPrimes, Real.log q)
        + 3 * (x.degOfDefinition : ℝ) * x.logDiffMell := by
  have hnn : ∀ q : ℕ, 0 ≤ Real.log q := by
    intro q
    rcases Nat.eq_zero_or_pos q with h | h
    · rw [h]; simp
    · exact Real.log_nonneg (by exact_mod_cast h)
  have hd : (x.degOfDefinition : ℝ) ≠ 0 := by
    have := x.degOfDefinition_pos
    positivity
  set D : ℝ := Real.log ((|NumberField.discr x.rep.toSSCurve.fld| : ℤ) : ℝ) with hD
  have hrhs : 3 * (x.degOfDefinition : ℝ) * x.logDiffMell = 3 * D := by
    rw [logDiffMell, ← hD]
    field_simp
  -- 分岐する素数の分
  have hdisc : (∑ q ∈ (NumberField.discr x.rep.toSSCurve.fld).natAbs.primeFactors, Real.log q)
      ≤ D := by
    have hpos : 0 < (NumberField.discr x.rep.toSSCurve.fld).natAbs :=
      Int.natAbs_pos.2 (NumberField.discr_ne_zero (K := x.rep.toSSCurve.fld))
    refine le_trans (sum_log_primeFactors_le _ hpos) (le_of_eq ?_)
    have hcast : ((|NumberField.discr x.rep.toSSCurve.fld| : ℤ) : ℝ)
        = (((NumberField.discr x.rep.toSSCurve.fld).natAbs : ℕ) : ℝ) := by
      rw [← Int.natCast_natAbs, Int.cast_natCast]
    rw [hD, hcast]
  -- `d` 以下の素数の分
  have htheta : (∑ q ∈ (Finset.range (x.degOfDefinition + 1)).filter Nat.Prime, Real.log q)
      ≤ 2 * D := by
    have hlog : (∑ q ∈ (Finset.range (x.degOfDefinition + 1)).filter Nat.Prime, Real.log q)
        = Real.log (primorial x.degOfDefinition) := by
      rw [primorial]
      push_cast
      rw [Real.log_prod]
      intro q hq
      exact_mod_cast (Finset.mem_filter.1 hq).2.ne_zero
    rw [hlog, hD]
    exact ABC3.Found.GaloisRep.log_primorial_le_two_log_discr x.rep.toSSCurve.fld
  have hu1 := sum_union_le_of_nonneg
    (x.rep.badPrimes ∪ (NumberField.discr x.rep.toSSCurve.fld).natAbs.primeFactors)
    ((Finset.range (x.degOfDefinition + 1)).filter Nat.Prime) (fun q : ℕ => Real.log q) hnn
  have hu2 := sum_union_le_of_nonneg x.rep.badPrimes
    ((NumberField.discr x.rep.toSSCurve.fld).natAbs.primeFactors) (fun q : ℕ => Real.log q) hnn
  rw [hrhs, ramPrimes]
  linarith

end RealizedClass

/-! ## ★★★★★★★★★★★★★★`quotLCyclic` 欄 -/

/-- ★★★★**Vélu の商を `SSCurve` として**——定義体は変えない。 -/
noncomputable def quotSSCurve (E : SSCurve) (S : Finset (E.fld × E.fld))
    (hell : (veluQuotientFull E.W S).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 E.fld), SemistableAt p (veluQuotientFull E.W S)) :
    SSCurve :=
  { K := E.K, W := veluQuotientFull E.W S, isEll := hell, ss := hss }

@[simp] theorem quotSSCurve_fld (E : SSCurve) (S : Finset (E.fld × E.fld)) (hell) (hss) :
    (quotSSCurve E S hell hss).fld = E.fld := rfl

/-- ★★★★★**`y` は `x` の `l`-巡回部分群による Vélu の商の類である**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`S` は `H∖{O}` の座標の集合であり、`|S| + 1 = l` が「位数 `l` の巡回群」の内容である。
☆`S` が実際に部分群の座標集合であることは、ここでは要求していない
——`degInf_quotLCyclic` 欄がそれを要求する側である。 -/
def IsQuotClassJ (x : RealizedClass) (l : ℕ) (y : ℂ) : Prop :=
  ∃ (S : Finset (x.rep.toSSCurve.fld × x.rep.toSSCurve.fld))
    (hell : (veluQuotientFull x.rep.toSSCurve.W S).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 x.rep.toSSCurve.fld),
      SemistableAt p (veluQuotientFull x.rep.toSSCurve.W S)),
    S.card + 1 = l ∧ (quotSSCurve x.rep.toSSCurve S hell hss).j = y

open scoped Classical in
/-- ★★★★★**`quotLCyclic` 欄**——商の類（無ければ自分自身）。 -/
noncomputable def quotLCyclicJ (x : RealizedClass) (l : ℕ) : RealizedClass :=
  if h : ∃ y : RealizedClass, IsQuotClassJ x l y.1 then h.choose else x

open scoped Classical in
theorem quotLCyclicJ_spec (x : RealizedClass) (l : ℕ)
    (h : ∃ y : RealizedClass, IsQuotClassJ x l y.1) :
    IsQuotClassJ x l (quotLCyclicJ x l).1 := by
  rw [quotLCyclicJ, dif_pos h]
  exact h.choose_spec

/-! ## ★出典の紐付け(`.src`)——★★条つき（半安定に制限した形） -/

def SSCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(EllModuliData の Curve 欄——半安定な楕円曲線の族)",
    sectionId := "genell-prop-3-4" }

def htFalt_congr_j.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(j が同じ SSCurve は ht^Falt も同じ。★無条件)",
    sectionId := "genell-prop-3-4" }

def faltingsHeightJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(EllModuliData の faltingsHeight 欄——j の函数として)",
    sectionId := "genell-prop-3-4" }

def SSCurve.localHeightAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(素点ごとの局所高さ——半安定なので v_p(Δ_min))",
    sectionId := "genell-lemma-3-7" }

def SSCurve.primeToLocalHeights_of_lt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(primeToLocalHeights_of_lt 欄——SSCurve の言葉で。★無条件)",
    sectionId := "genell-lemma-3-7" }

def torsionExt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(torsionExt 欄——半安定に制限したので恒等写像でよい)",
    sectionId := "genell-thm-3-8" }

def torsionExt.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆代償の記録: Curve := SSCurve(半安定なものだけ)と決めたので、" ++
       "原文 p.20 の L′(3･5 捧れを有理化する次数 23040 の拡大)の段は不要になるが、" ++
       "結論は「半安定な曲線について」に限定される") 3 ]

def quotSSCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商を SSCurve として——定義体は変えない)",
    sectionId := "genell-lemma-3-5" }

def quotLCyclicJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(quotLCyclic 欄——l-巡回部分群による商の類)",
    sectionId := "genell-lemma-3-5" }

def RealizedClass.logDiffMell.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(logDiffMell 欄——log|disc L| / d)",
    sectionId := "genell-cor-4-3" }

def RealizedClass.sum_log_ramPrimes_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(sum_log_ramPrimes_le 欄。★無条件)",
    sectionId := "genell-cor-4-3" }

def RealizedClass.sum_log_ramPrimes_le.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆原文は『by considering the trace of an extension of number fields』" ++
       "(分岐指数より 1 小さい重複度)で済ませているが、本形式化は" ++
       "Hermite-Minkowski(NumberField.abs_discr_ge)と Chebyshev(primorial_le_four_pow)で" ++
       "置き換えた。係数 3 の内訳は 1 単位(分岐する素数)＋2 単位(d 以下の素数)である") 4 ]

def RealizedClass.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(Curve 欄——実現される類と、次数最小の代表元)",
    sectionId := "genell-cor-4-3" }

def RealizedClass.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★測定: EllModuliData の logDiffMell は類の函数でなければならないが、" ++
       "log-diff は定義体の分岐で決まる量である。同じ j を持つ曲線でも定義体は" ++
       "いくらでも大きく取れるので、sum_log_ramPrimes_le 欄は「勝手な定義体の曲線」に" ++
       "対しては破れる。そこで Curve を実現される類とし、各類に定義体の次数が" ++
       "最小の代表元を取る——原文 Corollary 4.3 の仮説" ++
       "『L is a minimal field of definition of the point [E_L]』そのものである") 4 ]

def GaloisFiniteJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(GaloisFinite 欄——有限集合として取る)",
    sectionId := "genell-lemma-3-7" }

def GaloisFiniteJ.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆原文の Galois-finite(Example 1.3, (i))は「有限個の Galois 軌道の合併」であり、" ++
       "代数的数の Galois 軌道は有限なので集合として有限である。" ++
       "逆は成り立たないので、これは原文より弱い述語である" ++
       "——GaloisFinite は結論の側に現れるので、弱い方を取ると主張は強くなる") 3 ]

def noMultRedExcJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(noMultRedExc 欄——DegCurve なら空集合でよい)",
    sectionId := "genell-lemma-3-7" }

def MinimalFieldJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(MinimalField 欄——定義体の次数が最小)",
    sectionId := "genell-cor-4-3" }

def DegCurve.badPrimes.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(badPrimes 欄——S∘ のうち S 以外の部分)",
    sectionId := "genell-cor-4-3" }

def DegCurve.sum_log_badPrimes_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(sum_log_badPrimes_le 欄。★無条件)",
    sectionId := "genell-cor-4-3" }

def DegCurve.sum_localHt_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(sum_localHt_eq 欄——∑ h_j·log(p_j) = 23040·d·deg∞。★無条件)",
    sectionId := "genell-cor-4-3" }

def DegCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(Curve 欄——乗法還元を持つ半安定曲線)",
    sectionId := "genell-cor-4-3" }

def DegCurve.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★測定(2026-08-31、第 745): EllModuliData の multCard_pos・localHt_pos・" ++
       "multPrime_prime・sum_localHt_eq は、どの E : Curve に対しても deg∞(cls E) > 0 を" ++
       "強制する(Check/GenEll/EllModuliDegInfPos.lean)。" ++
       "したがって Curve 欄は「至る所良還元」の曲線を含めない") 3 ]

def northcottJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(northcott 欄。★無条件)",
    sectionId := "genell-prop-3-4" }

def degLeJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(EllModuliData の degLe 欄——M_ell(ℚ̄)^{≤d})",
    sectionId := "genell-prop-3-4" }

def degInfJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(EllModuliData の degInf 欄——j の函数として)",
    sectionId := "genell-prop-3-4" }

def degInfJ_le_faltingsHeightJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(類の水準で deg∞ ≤ 12·ht^Falt + A。★無条件)",
    sectionId := "genell-prop-3-4" }

def SSCurve.htFalt_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(SSCurve の上で ht^Falt は変数変換で不変。★無条件)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
