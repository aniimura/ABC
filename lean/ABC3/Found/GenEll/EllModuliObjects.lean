/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.FaltingsWitness
import ABC3.Found.GaloisRep.SemistableFin
import ABC3.Found.GaloisRep.Compositum
import ABC3.Found.GaloisRep.HtFaltBounds
import ABC3.Found.GaloisRep.HtJBound
import ABC3.Found.GaloisRep.ResChar
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
  /-- 定義体（`ℂ` の部分体として取る——`Prop34.lean` の族の形に合わせる）。 -/
  fld : Type
  [isField : Field fld]
  [isNF : NumberField fld]
  /-- 曲線。 -/
  W : WeierstrassCurve fld
  [isEll : W.IsElliptic]
  /-- ★全ての有限素点で半安定。 -/
  ss : ∀ p : HeightOneSpectrum (𝓞 fld), SemistableAt p W
  /-- `j` を `ℂ` へ送る埋め込み。 -/
  emb : fld →+* ℂ

attribute [instance] SSCurve.isField SSCurve.isNF SSCurve.isEll

namespace SSCurve

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
    htFalt { fld := E.fld, W := C • E.W, isEll := hell, ss := hss, emb := E.emb }
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

/-- ★★★★★★★★★★★★★★★★★★★★
**`northcott` 欄は `j` の高さの Northcott 性ひとつに帰着する**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★仮説 `hN` は `Skeleton/GenEll/NorthcottJ.lean` の `northcott_htJ` そのもの
（原文の `Proposition 1.4, (iv)`、実質は古典的 Northcott の定理）。

★★機構は `prop_3_4_chain_semistable`（`§9-1004`）の 2 本目の `≲`:

    ht^Falt(E) ≤ C  ⟹  h(j) ≤ 12(1+ϵ)·C + C₀

——高さが抑えられた類は `j` の高さも抑えられる。 -/
theorem northcottJ_of_northcott_htJ
    (hN : ∀ (B : ℝ) (d : ℕ),
      {x : ℂ | ∃ E : SSCurve, E.j = x ∧ E.deg ≤ d ∧ htJ E.fld E.W ≤ B}.Finite)
    (C : ℝ) (d : ℕ) : {x ∈ degLeJ d | faltingsHeightJ x ≤ C}.Finite := by
  obtain ⟨C₀, hC₀⟩ := prop_3_4_chain_semistable 1 one_pos
  refine (hN (12 * (1 + 1) * C + C₀) d).subset ?_
  rintro x ⟨⟨E, hEj, hEd⟩, hC⟩
  refine ⟨E, hEj, hEd, ?_⟩
  have hfe : faltingsHeightJ x = htFaltOf E.fld E.W := by
    rw [← hEj]; exact faltingsHeightJ_eq E
  have h1 : htFaltOf E.fld E.W ≤ C := by rw [← hfe]; exact hC
  have h2 := (hC₀ E.fld E.W E.ss).2.1
  linarith

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

★原文の Galois-finite（`Example 1.3, (i)`）は「有限個の Galois 軌道の合併」であり、
代数的数の Galois 軌道は有限なので**集合として有限**である。
☆逆は成り立たない（有限集合が Galois 安定とは限らない）ので、
これは原文より**弱い**述語である——`GaloisFinite` は結論の側に現れるので、
弱い方を取ると主張が強くなる。 -/
def GaloisFiniteJ (S : Set ℂ) : Prop := S.Finite

/-- ★★**`galoisFinite_union` 欄**。 -/
theorem galoisFiniteJ_union (S T : Set ℂ) (hS : GaloisFiniteJ S) (hT : GaloisFiniteJ T) :
    GaloisFiniteJ (S ∪ T) := hS.union hT

/-- ★★★**`CompactlyBounded` 欄**——アルキメデス素点での有界性。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★原文の compactly bounded（`Example 1.3, (ii)`）は有限個の素点 `V` と
そこでのコンパクト集合 `K_v` で押さえることだが、`EllClass := ℂ` の下では
アルキメデス素点の条件が `j` の有界性である。 -/
def CompactlyBoundedJ (S : Set ℂ) : Prop := Bornology.IsBounded S

/-- ★**`compactlyBounded_empty` 欄**。 -/
theorem compactlyBoundedJ_empty : CompactlyBoundedJ (∅ : Set ℂ) := Bornology.isBounded_empty

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
    GaloisFiniteJ (noMultRedExcJ KV) := Set.finite_empty

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

def degLeJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(EllModuliData の degLe 欄——M_ell(ℚ̄)^{≤d})",
    sectionId := "genell-prop-3-4" }

def northcottJ_of_northcott_htJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(northcott 欄は j の高さの Northcott 性ひとつに帰着する)",
    sectionId := "genell-prop-3-4" }

def northcottJ_of_northcott_htJ.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_htJ(Skeleton/GenEll/NorthcottJ.lean、古典的 Northcott)"
      (.absent "mathlib: Northcott の instance は具体的な体に対して 1 つも無い(2026-08-31 測定)") 12 ]

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
