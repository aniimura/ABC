/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Analysis.SpecialFunctions.Elliptic.Weierstrass
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Basic
import ABC3.Found.GenEll.LatticeCurve
import ABC3.Found.GenEll.WeierstrassODE
import ABC3.Found.GenEll.Velu
import ABC3.Found.GenEll.PointVariableChange
import ABC3.Meta.Claim
import ABC3.Found.GenEll.Uniformization.Phi

/-!
# 一様化 —— 部分群の原像・群同型 `ℂ/Λ ≅ E(ℂ)`・階数 1 の部分群・指数

☆`Found/GenEll/Uniformization.lean`(292 KB / 325 宣言)を
**ファイル内の見出し**で割った 1 枚である(2026-09-03、第 1456)。
★論文のセクションでは割れない——このファイルは [GenEll] §3 の
`Lemma 3.5` と `Proposition 3.4` の 2 項目しか持たず、割っても 146 KB のままだからである。
☆読む順序は `Basic → VeluAnalytic → Surjective → AdditionEntry → AdditionODE
→ FilledPole → AdditionFormula → Phi → GroupIso → Sublattice → G2G3 → Assemble`。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★部分群の原像 -/

open scoped Classical in
/-- ★★★★★★`Φ(0) = 0`。 -/
@[simp] theorem uniformMap_zero (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) :
    uniformMap P hΔ 0 = 0 :=
  uniformMap_of_mem P hΔ P.lattice.zero_mem

open scoped Classical in
/-- ★★★★★★★★`Φ(−z) = −Φ(z)`。 -/
theorem uniformMap_neg (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ) :
    uniformMap P hΔ (-z) = -uniformMap P hΔ z := by
  have h := uniformMap_add P hΔ z (-z)
  rw [add_neg_cancel, uniformMap_zero] at h
  exact (neg_eq_of_add_eq_zero_right h.symm).symm

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★**部分群 `H ⊆ E(ℂ)` の原像は ℂ の部分群**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★`Φ` が群準同型（第 661）だから。☆これが `Lemma 3.5` の `Λ′` である。 -/
noncomputable def preimageSubgroup (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) : AddSubgroup ℂ where
  carrier := {z : ℂ | uniformMap P hΔ z ∈ H}
  add_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq, uniformMap_add] at *
    exact H.add_mem ha hb
  zero_mem' := by
    simp only [Set.mem_setOf_eq, uniformMap_zero]
    exact H.zero_mem
  neg_mem' := by
    intro a ha
    simp only [Set.mem_setOf_eq, uniformMap_neg] at *
    exact H.neg_mem ha

open scoped Classical in
@[simp] theorem mem_preimageSubgroup (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) (z : ℂ) :
    z ∈ preimageSubgroup P hΔ H ↔ uniformMap P hΔ z ∈ H := Iff.rfl

open scoped Classical in
/-- ★★★★★★★★★★**`Λ ⊆ Λ′`**——格子は `Φ` で `0` に落ちるから。 -/
theorem lattice_le_preimageSubgroup (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) {z : ℂ} (hz : z ∈ P.lattice) :
    z ∈ preimageSubgroup P hΔ H := by
  simp only [mem_preimageSubgroup, uniformMap_of_mem P hΔ hz]
  exact H.zero_mem

open scoped Classical in
/-- ★★★★★★★★★★★★★★**`Φ` は `Λ` を法として単射**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★第 624 の単射性を `Φ` の言葉に直したもの。☆`Λ′/Λ → H` の単射性に要る。 -/
theorem sub_mem_lattice_of_uniformMap_eq (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {z w : ℂ} (h : uniformMap P hΔ z = uniformMap P hΔ w) : z - w ∈ P.lattice := by
  by_cases hz : z ∈ P.lattice
  · by_cases hw : w ∈ P.lattice
    · exact P.lattice.sub_mem hz hw
    · rw [uniformMap_of_mem P hΔ hz, uniformMap_of_notMem P hΔ hw] at h
      exact absurd h.symm (by simp)
  · by_cases hw : w ∈ P.lattice
    · rw [uniformMap_of_notMem P hΔ hz, uniformMap_of_mem P hΔ hw] at h
      exact absurd h (by simp)
    · rw [uniformMap_of_notMem P hΔ hz, uniformMap_of_notMem P hΔ hw] at h
      have hx : latticePointX P z = latticePointX P w := by
        injection h with hx hy
      have hy : latticePointY P z = latticePointY P w := by
        injection h with hx hy
      have hpx : P.weierstrassP z = P.weierstrassP w := hx
      have hpy : P.derivWeierstrassP z = P.derivWeierstrassP w := by
        simp only [latticePointY] at hy
        linear_combination 2 * hy
      refine mem_lattice_of_shift_eq P (z - w) hw ?_ ?_ ?_
      · rw [show w + (z - w) = z by ring]; exact hz
      · rw [show w + (z - w) = z by ring]; exact hpx
      · rw [show w + (z - w) = z by ring]; exact hpy

def preimageSubgroup.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(部分群 H ⊆ E(ℂ) の原像は ℂ の部分群——これが Λ′ である)",
    sectionId := "genell-lemma-3-5" }

def sub_mem_lattice_of_uniformMap_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Φ は Λ を法として単射。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★群同型 `ℂ/Λ ≅ E(ℂ)` -/

open scoped Classical in
/-- ★★★★★★★★★★`Φ(z) = 0 ⟺ z ∈ Λ`——`Φ` の核はちょうど `Λ`。 -/
theorem uniformMap_eq_zero_iff (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ) :
    uniformMap P hΔ z = 0 ↔ z ∈ P.lattice := by
  refine ⟨fun h => ?_, fun h => uniformMap_of_mem P hΔ h⟩
  by_contra hz
  rw [uniformMap_of_notMem P hΔ hz] at h
  exact absurd h (by simp)

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**一様化写像を加法群準同型として束ねたもの**。

★中身は第 661 の `uniformMap_add`。 -/
noncomputable def uniformHom (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) :
    ℂ →+ (latticeCurve P).toAffine.Point :=
  AddMonoidHom.mk' (uniformMap P hΔ) (uniformMap_add P hΔ)

open scoped Classical in
@[simp] theorem uniformHom_apply (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ) :
    uniformHom P hΔ z = uniformMap P hΔ z := rfl

open scoped Classical in
/-- ★★★★★★★★★★★★`Φ` の核は `Λ`。 -/
theorem ker_uniformHom (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) :
    (uniformHom P hΔ).ker = P.lattice.toAddSubgroup := by
  ext z
  simp only [AddMonoidHom.mem_ker, uniformHom_apply, uniformMap_eq_zero_iff,
    Submodule.mem_toAddSubgroup]

open scoped Classical in
/-- ★★★★★★★★★★★★★★**`Φ` は全射**——第 604 の `latticePoint_surjective` を
`Point` の言葉に直したもの。 -/
theorem uniformMap_surjective (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) :
    Function.Surjective (uniformMap P hΔ) := by
  intro Q
  cases Q with
  | zero => exact ⟨0, uniformMap_zero P hΔ⟩
  | some x y h =>
      obtain ⟨z, hz, hx, hy⟩ := latticePoint_surjective P x y h.left
      refine ⟨z, ?_⟩
      rw [uniformMap_of_notMem P hΔ hz]
      exact point_some_congr hx hy

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**一様化定理——`ℂ/Λ ≅ E(ℂ)` は加法群の同型**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★三つの部品がすべて揃った:

* **全射** —— 第 603（`weierstrassP_surjective`）＋ 第 604
* **単射** —— 第 624（`mem_lattice_of_shift_eq`）＋ 第 662
* **準同型** —— 第 661（`uniformMap_add`）

★★★★★☆**どの部品も mathlib に無い**（`§9-1039` で測った通り）。
☆これで `Lemma 3.5` の「`l`-捻れの部分群 ↔ `Λ` を含む指数 `l` の格子」が
純粋に群論の言葉で書けるようになる。 -/
noncomputable def uniformEquiv (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) :
    (ℂ ⧸ P.lattice.toAddSubgroup) ≃+ (latticeCurve P).toAffine.Point :=
  AddEquiv.ofBijective
    (QuotientAddGroup.lift P.lattice.toAddSubgroup (uniformHom P hΔ)
      (fun z hz => uniformMap_of_mem P hΔ hz))
    ⟨by
      intro a b hab
      obtain ⟨z, rfl⟩ := QuotientAddGroup.mk_surjective a
      obtain ⟨w, rfl⟩ := QuotientAddGroup.mk_surjective b
      simp only [QuotientAddGroup.lift_mk, uniformHom_apply] at hab
      rw [QuotientAddGroup.eq, Submodule.mem_toAddSubgroup,
        show -z + w = -(z - w) by ring]
      exact P.lattice.neg_mem (sub_mem_lattice_of_uniformMap_eq P hΔ hab),
     by
      intro Q
      obtain ⟨z, hz⟩ := uniformMap_surjective P hΔ Q
      exact ⟨(z : ℂ ⧸ P.lattice.toAddSubgroup), by simpa using hz⟩⟩

open scoped Classical in
@[simp] theorem uniformEquiv_mk (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ) :
    uniformEquiv P hΔ (z : ℂ ⧸ P.lattice.toAddSubgroup) = uniformMap P hΔ z := rfl

def uniformEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(一様化定理——ℂ/Λ ≅ E(ℂ) は加法群の同型。★無条件)",
    sectionId := "genell-lemma-3-5" }

def uniformMap_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Φ は全射。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★階数 1 の部分群（`Lemma 3.5`） -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**位数 `l` の点 `Q` を生む `z₀`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`Φ` は全射（第 663）だから `Φ(z₀) = Q` なる `z₀` が取れ、`Q ≠ 0` だから
`z₀ ∉ Λ`、`l·Q = 0` だから `Φ(l z₀) = 0`、すなわち `l z₀ ∈ Λ`。
☆この `z₀` が「`Λ` に `1/l` の分母で足す元」である。 -/
theorem exists_generator_of_torsion (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} (hQ : Q ≠ 0) {l : ℕ} (hl : l • Q = 0) :
    ∃ z₀ : ℂ, z₀ ∉ P.lattice ∧ uniformMap P hΔ z₀ = Q ∧ (l : ℂ) * z₀ ∈ P.lattice := by
  obtain ⟨z₀, hz₀⟩ := uniformMap_surjective P hΔ Q
  refine ⟨z₀, fun hc => hQ ?_, hz₀, ?_⟩
  · rw [← hz₀, uniformMap_of_mem P hΔ hc]
  · have hzero : uniformMap P hΔ ((l : ℂ) * z₀) = 0 := by
      have h1 : ((l : ℂ) * z₀) = l • z₀ := by simp [nsmul_eq_mul]
      rw [h1, ← uniformHom_apply, map_nsmul, uniformHom_apply, hz₀, hl]
    exact (uniformMap_eq_zero_iff P hΔ _).1 hzero

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★
**巡回部分群の原像は `Λ′ = Λ + ℤ z₀`**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★`⊇` は準同型性（第 661）から、`⊆` は単射性（第 662）から:
`Φ(z) = k·Φ(z₀) = Φ(k z₀)` なら `z − k z₀ ∈ Λ`。

★★★★☆**これが `Lemma 3.5` の「階数 1」の内容である**——`E(ℂ)` の巡回部分群
`⟨Q⟩` は `Λ` に元を 1 つ添加した格子に対応する。 -/
theorem preimageSubgroup_zmultiples (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z₀ : ℂ) :
    preimageSubgroup P hΔ (AddSubgroup.zmultiples (uniformMap P hΔ z₀))
      = P.lattice.toAddSubgroup ⊔ AddSubgroup.zmultiples z₀ := by
  ext z
  simp only [mem_preimageSubgroup, AddSubgroup.mem_zmultiples_iff, AddSubgroup.mem_sup,
    Submodule.mem_toAddSubgroup]
  constructor
  · rintro ⟨k, hk⟩
    have hzk : uniformMap P hΔ (k • z₀) = uniformMap P hΔ z := by
      rw [← uniformHom_apply, map_zsmul, uniformHom_apply, hk]
    exact ⟨z - k • z₀, sub_mem_lattice_of_uniformMap_eq P hΔ hzk.symm, k • z₀,
      ⟨k, rfl⟩, by abel⟩
  · rintro ⟨y, hy, w, ⟨k, rfl⟩, rfl⟩
    refine ⟨k, ?_⟩
    have h2 : uniformMap P hΔ (y + k • z₀) = k • uniformMap P hΔ z₀ := by
      rw [← uniformHom_apply, map_add, map_zsmul, uniformHom_apply, uniformHom_apply,
        uniformMap_of_mem P hΔ hy, zero_add]
    exact h2.symm

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★
**`Λ ⊆ Λ′` かつ `l·Λ′ ⊆ Λ`**——`Λ′` が `Λ` の「`1/l` 倍の中」にある。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆これで `Λ′` は `Λ` と `(1/l)Λ` に挟まれ、有限指数の格子であることが確定する
（`Λ ⊆ Λ′ ⊆ (1/l)Λ`、`[(1/l)Λ : Λ] = l²`）。 -/
theorem smul_preimageSubgroup_le (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {H : AddSubgroup (latticeCurve P).toAffine.Point} {l : ℕ}
    (hH : ∀ Q ∈ H, l • Q = 0) {z : ℂ} (hz : z ∈ preimageSubgroup P hΔ H) :
    (l : ℂ) * z ∈ P.lattice := by
  have hzero : uniformMap P hΔ ((l : ℂ) * z) = 0 := by
    have h1 : ((l : ℂ) * z) = l • z := by simp [nsmul_eq_mul]
    rw [h1, ← uniformHom_apply, map_nsmul, uniformHom_apply]
    exact hH _ ((mem_preimageSubgroup P hΔ H z).1 hz)
  exact (uniformMap_eq_zero_iff P hΔ _).1 hzero

def exists_generator_of_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(位数 l の点 Q を生む z₀ の存在。★無条件)",
    sectionId := "genell-lemma-3-5" }

def preimageSubgroup_zmultiples.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(巡回部分群の原像は Λ′ = Λ + ℤz₀——階数 1 の内容。★無条件)",
    sectionId := "genell-lemma-3-5" }

def smul_preimageSubgroup_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(l·Λ′ ⊆ Λ——Λ′ は Λ と (1/l)Λ に挟まれる。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★指数 `[Λ′ : Λ] = |H|` -/

open scoped Classical in
/-- ★★★★★★★★★★`Φ` を `Λ′ → H` に制限したもの。 -/
noncomputable def preimageToH (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) :
    preimageSubgroup P hΔ H →+ H where
  toFun z := ⟨uniformMap P hΔ z.1, (mem_preimageSubgroup P hΔ H z.1).1 z.2⟩
  map_zero' := Subtype.ext (by simpa using uniformMap_zero P hΔ)
  map_add' := fun x y => Subtype.ext (by simpa using uniformMap_add P hΔ x.1 y.1)

open scoped Classical in
@[simp] theorem preimageToH_coe (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) (z : preimageSubgroup P hΔ H) :
    (preimageToH P hΔ H z : (latticeCurve P).toAffine.Point) = uniformMap P hΔ z.1 := rfl

open scoped Classical in
/-- ★★★★★★★★★★★★`Λ′ → H` は全射——`Φ` が全射（第 663）だから。 -/
theorem preimageToH_surjective (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) :
    Function.Surjective (preimageToH P hΔ H) := by
  rintro ⟨Q, hQ⟩
  obtain ⟨z, hz⟩ := uniformMap_surjective P hΔ Q
  refine ⟨⟨z, ?_⟩, Subtype.ext ?_⟩
  · rw [mem_preimageSubgroup, hz]; exact hQ
  · exact hz

open scoped Classical in
/-- ★★★★★★★★★★★★`Λ′ → H` の核は `Λ`（`Λ′` の中で見たもの）。 -/
theorem ker_preimageToH (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) :
    (preimageToH P hΔ H).ker
      = P.lattice.toAddSubgroup.addSubgroupOf (preimageSubgroup P hΔ H) := by
  ext z
  simp only [AddMonoidHom.mem_ker, AddSubgroup.mem_addSubgroupOf, Submodule.mem_toAddSubgroup]
  constructor
  · intro h
    exact (uniformMap_eq_zero_iff P hΔ _).1 (congrArg Subtype.val h)
  · intro h
    exact Subtype.ext (by simpa using uniformMap_of_mem P hΔ h)

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Λ′/Λ ≅ H`**——原像の商はもとの部分群と同型。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★全射（第 663）＋核 `= Λ`（本ブロック）＋第一同型定理。 -/
noncomputable def preimageQuotEquiv (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) :
    (preimageSubgroup P hΔ H
        ⧸ P.lattice.toAddSubgroup.addSubgroupOf (preimageSubgroup P hΔ H)) ≃+ H :=
  (QuotientAddGroup.quotientAddEquivOfEq (ker_preimageToH P hΔ H).symm).trans
    (QuotientAddGroup.quotientKerEquivOfSurjective (preimageToH P hΔ H)
      (preimageToH_surjective P hΔ H))

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`[Λ′ : Λ] = |H|`**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★☆**これで `Lemma 3.5` の「位数 `l` の部分群 ↔ 指数 `l` の格子」の
指数の部分が塞がった**——`H` が位数 `l` なら `Λ ⊆ Λ′` は指数 `l`。
☆残るのは `Λ′` の基底を取り出して行列式 `= l` を書くこと
（`Λ ⊆ Λ′ ⊆ (1/l)Λ` は第 664 で押さえてある）。 -/
theorem relIndex_preimageSubgroup (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) :
    P.lattice.toAddSubgroup.relIndex (preimageSubgroup P hΔ H) = Nat.card H := by
  rw [AddSubgroup.relIndex, AddSubgroup.index]
  exact Nat.card_congr (preimageQuotEquiv P hΔ H).toEquiv

def preimageQuotEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Λ′/Λ ≅ H——原像の商はもとの部分群と同型。★無条件)",
    sectionId := "genell-lemma-3-5" }

def relIndex_preimageSubgroup.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5([Λ′ : Λ] = |H|——指数は部分群の位数。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
