import ABC3.Found.FrdI.Prop21
import ABC3.Check.FrdI.TwistedFrobenioid

/-!
# ★★★[FrdI] `Proposition 2.1, (iii)` を「`d` 固定」で読むと**偽**である

原文 (FrdI p.44):
> (iii) C is of perfect type if and only if Ψ is an equivalence of categories.

★原文は命題の冒頭で `d ∈ ℕ≥1` を固定するが、`perfect object`
(`Definition 1.2, (iv)`、原文 p.23)は

原文 (FrdI p.23):
> C such that for every n ∈N≥1, it holds that every B ∈Ob(C) base-isomorphic to

と★**すべての `n ∈ ℕ≥1`** を量化する。★したがって `⟸` は 1 つの `d` からは出ない。

★★**ここではそれを「機械検証された事実」にする** ——
`(naiveFrob P F d).IsEquivalence` が成り立つのに `IsOfPerfectType P` が成り立たない
`(P, d)` を**実際に構成する**。

## ★構成

★★**`d = 1` で足りる。** `Definition 1.3, (ii)` が与える次数 1 の Frobenius 型射は
`Proposition 1.4, (iii)`(LB-invertible な pre-step は同型)により**同型**なので、
★**`Ψ_1` はつねに恒等関手と自然同型**であり、したがってつねに圏同値である
(`naiveFrob_one_isEquivalence`、★これは任意の Frobenioid についての一般の事実)。

★一方 `𝔽_ℕ`(`𝒟 = Discrete PUnit`、`Φ = ℕ`)は **perfect 型でない** ——
`n = 2` で perfect の第 2 条件が破れる(`ℕ` で `2a = 1` が解けない)。

★★★**したがって「`d` 固定」の読みは、`𝔽_ℕ` という原文自身の基本例で破れる。**
★しかも `Ψ_1 ≅ 𝟭` は**すべての** Frobenioid で成り立つので、
★★**「`d` 固定」の読みは perfect 型でない Frobenioid すべてで破れる。**

★★**さらに `d ≥ 2` の例も作る**(第 4 段) —— 「`d = 1` は退化だ」を塞ぐため。
`Φ := ℤ/2`(定数)では `3 • x = x`(全単射)・`2 • x = 0`(非全射)なので、
★**`Ψ_3` は圏同値だが `n = 2` で perfect が破れる**
(`prop_2_1_iii_fixed_degree_false_deg_three`)。

★これは `Found/FrdI/Prop21.lean` の `prop_2_1_iii` が
`∀ n : ℕ+, (naiveFrob P F n).IsEquivalence` という形を取っている理由である。

## ★機械検証の要点(`#print axioms` で確認済み、2026-08-17)

`prop_2_1_iii_fixed_degree_false` / `_general` / `_deg_three` はいずれも
`[propext, Classical.choice, Quot.sound]` のみに依存する(★`sorryAx` は無い)。
-/

namespace ABC3.Check.FrdI

open CategoryTheory ABC3.Found.FrdI

universe v u w u2 v2

/-! ## ★第 1 段 —— `Ψ_1` はつねに圏同値(一般) -/

section General

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-- ★**次数 1 の Frobenius 型射は同型** ——
Frobenius 型は LB-invertible、次数 1 かつ底同型なので pre-step であり、
`Proposition 1.4, (iii)` が「LB-invertible な pre-step は同型」を与える。 -/
theorem nfHom_one_isIso (F : FrobenioidCore P) (A : C) : IsIso (nfHom P F 1 A) :=
  prop_1_4_iii P F (nfHom P F 1 A) (nfHom_frobType P F 1 A).1
    ⟨nfHom_degFr P F 1 A, (nfHom_frobType P F 1 A).2⟩

/-- ★★**`Ψ_1` は恒等関手と自然同型** —— 単位 `α : 𝟭 𝒞 ⟶ Ψ_1` の各成分が同型。 -/
noncomputable def naiveFrobOneIso (F : FrobenioidCore P) : 𝟭 C ≅ naiveFrob P F 1 :=
  haveI : ∀ A : C, IsIso ((naiveFrobUnit P F 1).app A) := nfHom_one_isIso P F
  haveI : IsIso (naiveFrobUnit P F 1) := NatIso.isIso_of_isIso_app _
  asIso (naiveFrobUnit P F 1)

/-- ★★★**次数 1 の naive Frobenius 関手はつねに圏同値**。

★★**perfect 型かどうかに一切依らない** —— ここが「`d` 固定」の読みが破れる急所である。 -/
theorem naiveFrob_one_isEquivalence (F : FrobenioidCore P) :
    (naiveFrob P F 1).IsEquivalence :=
  Functor.isEquivalence_of_iso (naiveFrobOneIso P F)

end General

/-! ## ★第 2 段 —— `𝔽_ℕ` は perfect 型でない

★`𝒟 = Discrete PUnit`、`Φ = ℕ`(定数)の**捻れていない**素の `𝔽_ℕ` を使う
(`Check/FrdI/TwistedFrobenioid.lean` の `efP` / `efA`)。

★`ElemFrobCat.div_comp` より `(a,f) ≫ (m,d)` の `Div` は `m + d·a` なので、
次数 2 の Frobenius 型射 `(0,2)` を右から掛けると `Div` は**偶数**しか作れない。
★したがって `Div = 1` の pre-step `(1,1)` を「割る」ことができない。 -/

/-- ★`𝔽_ℕ` の `FrobenioidCore`(`Proposition 1.5`)。 -/
theorem efF : FrobenioidCore efP :=
  elemFrob_frobenioidCore Cx2Phi cx2_totEpi (fun _ => isPreDivisorial_nat)

/-- ★次数 2 の Frobenius 型射 `(Div 0, deg 2)`。 -/
def efFrob2 : efA ⟶ efA := ⟨𝟙 _, (0 : ℕ), 2⟩

/-- ★`Div = 1` の pre-step。 -/
def efStep1 : efA ⟶ efA := ⟨𝟙 _, (1 : ℕ), 1⟩

theorem efFrob2_frobType : IsFrobeniusType efP efFrob2 := by
  refine (elemFrob_frobeniusType_iff Cx2Phi cx2_totEpi (fun _ => isPreDivisorial_nat)
    efFrob2).mpr ⟨?_, ?_⟩
  · show toChar (ElemFrobCat.Hom.div efFrob2) = 0
    show toChar (0 : ℕ) = 0
    exact map_zero _
  · show IsIso (ElemFrobCat.Hom.base efFrob2)
    show IsIso (𝟙 _)
    infer_instance

theorem efFrob2_deg : efP.degFr efFrob2 = 2 := rfl

theorem efStep1_preStep : IsPreStep efP efStep1 := by
  refine ⟨rfl, ?_⟩
  show IsIso (ElemFrobCat.Hom.base efStep1)
  show IsIso (𝟙 _)
  infer_instance

/-- ★★★**`𝔽_ℕ` は perfect 型でない** —— `n = 2` で perfect の第 2 条件が破れる。

★`efFrob2 ≫ efStep1` の `Div` は `1`、`ψ ≫ efFrob2` の `Div` は `2 · ψ.div` で、
★**`ℕ` では `2a = 1` に解が無い**。 -/
theorem ef_not_perfectType : ¬ IsOfPerfectType efP := by
  intro h
  obtain ⟨-, huniq⟩ := h efA 2
  obtain ⟨ψ, ⟨-, hsq⟩, -⟩ := huniq efA efA efA efA efFrob2 efFrob2
    efFrob2_frobType efFrob2_deg efFrob2_frobType efFrob2_deg
    ⟨Iso.refl _⟩ ⟨Iso.refl _⟩ efStep1 efStep1_preStep
  have hd := congrArg ElemFrobCat.Hom.div hsq
  rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp] at hd
  simp only [efFrob2, efStep1, constPhi_map, zero_add, smul_zero, add_zero] at hd
  rw [show ((2 : ℕ+) : ℕ) = 2 from rfl, two_nsmul] at hd
  obtain ⟨a, ha⟩ : ∃ a : ℕ, (1 : ℕ) + ((1 : ℕ+) : ℕ) • (0 : ℕ) = 0 + (a + a) :=
    ⟨ElemFrobCat.Hom.div ψ, hd⟩
  simp only [smul_zero, add_zero, zero_add] at ha
  omega

/-! ## ★第 3 段 —— 結論 -/

/-- ★★★**[FrdI] `Proposition 2.1, (iii)` の「`d` 固定」の読みは偽である**。

★**`d` を固定したまま `(naiveFrob P F d).IsEquivalence → IsOfPerfectType P` と読むと、
`d = 1`・`P = 𝔽_ℕ` で反例になる。**

★★これが `Found/FrdI/Prop21.lean` の `prop_2_1_iii` を
`IsOfPerfectType P ↔ ∀ n : ℕ+, (naiveFrob P F n).IsEquivalence`
という形で述べた理由である。 -/
theorem prop_2_1_iii_fixed_degree_false :
    ∃ d : ℕ+, (naiveFrob efP efF d).IsEquivalence ∧ ¬ IsOfPerfectType efP :=
  ⟨1, naiveFrob_one_isEquivalence efP efF, ef_not_perfectType⟩

/-- ★**より強い形** —— 「`d` 固定」の読みは、perfect 型でない Frobenioid
**すべて**で破れる(`Ψ_1 ≅ 𝟭` はつねに成り立つから)。 -/
theorem prop_2_1_iii_fixed_degree_false_general
    {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
    {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
    (hnp : ¬ IsOfPerfectType P) :
    ∃ d : ℕ+, (naiveFrob P F d).IsEquivalence ∧ ¬ IsOfPerfectType P :=
  ⟨1, naiveFrob_one_isEquivalence P F, hnp⟩


/-! ## ★第 4 段 —— `d ≥ 2` の例(`d = 3`、`Φ = ℤ/2`)

★★**「`d = 1` は退化だ」という反論を塞ぐ。**

★`Φ := ℤ/2`(定数)を取ると `3 • x = x`(**全単射**)・`2 • x = 0`(**非全射**)なので、
★**`Ψ_3` は圏同値だが `n = 2` で perfect の第 2 条件が破れる。**

★`ℤ/2` は群なので sharp ではないが、`Definition 1.1, (iv)` が要求するのは
**`Φ` が pre-divisorial** であることで(除子モノイドは `Φ^char`)、
`Proposition 1.5` はそこから `𝔽_Φ` が Frobenioid であることを与える。
★`ℤ` について同じ形の実例が既にある(`isPreDivisorial_int`)。 -/

theorem range_toGp_zmod2 : Set.range (toGp (ZMod 2)) = Set.univ := by
  ext x
  simp only [Set.mem_univ, iff_true]
  induction x using AddLocalization.ind with
  | _ p =>
    obtain ⟨a, b⟩ := p
    refine ⟨a - (b : ZMod 2), ?_⟩
    rw [toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists]
    exact ⟨0, by simp⟩

theorem isOfCharacteristicType_zmod2 : IsOfCharacteristicType (ZMod 2) := by
  intro a b _
  refine ⟨b - a, ⟨⟨⟨b - a, a - b, by ring, by ring⟩, rfl⟩, by ring⟩, ?_⟩
  rintro y ⟨-, hy⟩
  rw [hy]; ring

theorem isPreDivisorial_zmod2 : IsPreDivisorial (ZMod 2) :=
  ⟨isIntegralMonoid_of_isCancelAdd (ZMod 2),
   isSaturatedMonoid_of_range_eq_univ (ZMod 2) range_toGp_zmod2,
   isOfCharacteristicType_zmod2⟩

/-- ★★`ℤ/2` では **`3` 倍は恒等**(したがって全単射)。 -/
theorem zmod2_three_nsmul : ∀ m : ZMod 2, (3 : ℕ) • m = m := by decide

/-- ★★`ℤ/2` では **`2` 倍は `0`**(したがって非全射 —— `1` は像に無い)。 -/
theorem zmod2_two_nsmul : ∀ m : ZMod 2, (2 : ℕ) • m = 0 := by decide

theorem zmod2_add_self : ∀ m : ZMod 2, m + m = 0 := by decide

/-- ★一対象の底圏上の定数モノイド `ℤ/2`。 -/
abbrev Z2Phi : MonoidOn.{0, 0, 0} (Discrete PUnit) := constPhi (ZMod 2)

/-- ★`𝔽_{ℤ/2}`。 -/
abbrev z2P : PreFrobenioid (ElemFrobCat Z2Phi) Z2Phi.charOn :=
  elemPreFrobenioid Z2Phi cx2_totEpi (fun _ => isPreDivisorial_zmod2)

theorem z2F : FrobenioidCore z2P :=
  elemFrob_frobenioidCore Z2Phi cx2_totEpi (fun _ => isPreDivisorial_zmod2)

/-- ★`Φ` が定数なので `Φ(X)` は `ℤ/2` そのもの。★消約の instance を通す。 -/
instance z2Val_isCancelAdd (X : Discrete PUnit) : IsCancelAdd (Z2Phi.val X) :=
  inferInstanceAs (IsCancelAdd (ZMod 2))

/-- ★`Div` を `ℤ/2` の元として取り出す(型を揃えるため)。 -/
theorem z2_div_exists {X Y : ElemFrobCat Z2Phi} (f : X ⟶ Y) :
    ∃ m : ZMod 2, ElemFrobCat.Hom.div f = m := ⟨_, rfl⟩

/-- ★底圏が一点なので、任意の 2 対象の間に射を作れる。 -/
def z2Hom (X Y : ElemFrobCat Z2Phi) (m : ZMod 2) (e : ℕ+) : X ⟶ Y :=
  ⟨eqToHom (Subsingleton.elim _ _), m, e⟩

@[simp] theorem z2Hom_div (X Y : ElemFrobCat Z2Phi) (m : ZMod 2) (e : ℕ+) :
    ElemFrobCat.Hom.div (z2Hom X Y m e) = m := rfl

@[simp] theorem z2Hom_deg (X Y : ElemFrobCat Z2Phi) (m : ZMod 2) (e : ℕ+) :
    ElemFrobCat.Hom.deg (z2Hom X Y m e) = e := rfl

/-- ★`nfHom` の次数は `3`(`Hom.deg` の言葉で)。 -/
theorem z2_nfHom_deg (A : ElemFrobCat Z2Phi) :
    ElemFrobCat.Hom.deg (nfHom z2P z2F 3 A) = 3 :=
  nfHom_degFr z2P z2F 3 A

/-- ★`ℤ/2` では**奇数倍は恒等**。 -/
theorem zmod2_odd_nsmul (k : ℕ) (hk : k % 2 = 1) (m : ZMod 2) : k • m = m := by
  rw [nsmul_eq_mul, ← ZMod.natCast_mod, hk]
  simp

/-- ★`ℤ/2` では**偶数倍は `0`**。 -/
theorem zmod2_even_nsmul (k : ℕ) (hk : k % 2 = 0) (m : ZMod 2) : k • m = 0 := by
  rw [nsmul_eq_mul, ← ZMod.natCast_mod, hk]
  simp

/-- ★★**次数 3 の `nfHom` は mono** —— `ℤ/2` で `3` 倍が単射だから。 -/
theorem z2_nfHom_mono (A : ElemFrobCat Z2Phi) : Mono (nfHom z2P z2F 3 A) := by
  refine ⟨fun {X} f g hfg => ?_⟩
  have hdiv := congrArg ElemFrobCat.Hom.div hfg
  have hdeg := congrArg ElemFrobCat.Hom.deg hfg
  rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp] at hdiv
  rw [ElemFrobCat.degFr_comp, ElemFrobCat.degFr_comp, z2_nfHom_deg] at hdeg
  obtain ⟨fd, hfd⟩ := z2_div_exists f
  obtain ⟨gd, hgd⟩ := z2_div_exists g
  obtain ⟨c, hc⟩ := z2_div_exists (nfHom z2P z2F 3 A)
  simp only [constPhi_map, hfd, hgd, hc] at hdiv
  have hodd : ((ElemFrobCat.Hom.deg (nfHom z2P z2F 3 A) : ℕ+) : ℕ) % 2 = 1 := by
    rw [z2_nfHom_deg]; rfl
  have hdiv2 : c + ((ElemFrobCat.Hom.deg (nfHom z2P z2F 3 A) : ℕ+) : ℕ) • fd
      = c + ((ElemFrobCat.Hom.deg (nfHom z2P z2F 3 A) : ℕ+) : ℕ) • gd := hdiv
  rw [zmod2_odd_nsmul _ hodd, zmod2_odd_nsmul _ hodd] at hdiv2
  exact ElemFrobCat.Hom.ext (Subsingleton.elim _ _)
    (by rw [hfd, hgd, add_left_cancel hdiv2]) (mul_left_cancel hdeg)

theorem z2_naiveFrob_faithful : (naiveFrob z2P z2F 3).Faithful := by
  constructor
  intro A B f g hfg
  have hfg' : nfMap z2P z2F 3 f = nfMap z2P z2F 3 g := hfg
  haveI := z2_nfHom_mono B
  refine (cancel_mono (nfHom z2P z2F 3 B)).mp ?_
  rw [nfMap_sq, nfMap_sq, hfg']

theorem z2_naiveFrob_full : (naiveFrob z2P z2F 3).Full := by
  constructor
  intro A B g
  obtain ⟨gd, hgd⟩ := z2_div_exists g
  obtain ⟨c₁, hc₁⟩ := z2_div_exists (nfHom z2P z2F 3 A)
  obtain ⟨c₂, hc₂⟩ := z2_div_exists (nfHom z2P z2F 3 B)
  refine ⟨z2Hom A B (gd + ((ElemFrobCat.Hom.deg g : ℕ) • c₁) + c₂)
    (ElemFrobCat.Hom.deg g), ?_⟩
  show nfMap z2P z2F 3 _ = g
  refine (nfMap_uniq z2P z2F 3 _ g ?_).symm
  refine ElemFrobCat.Hom.ext (Subsingleton.elim _ _) ?_ ?_
  · rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp]
    simp only [constPhi_map, z2Hom_div, hgd, hc₁, hc₂]
    have hodd : ((ElemFrobCat.Hom.deg (nfHom z2P z2F 3 B) : ℕ+) : ℕ) % 2 = 1 := by
      rw [z2_nfHom_deg]; rfl
    show c₂ + ((ElemFrobCat.Hom.deg (nfHom z2P z2F 3 B) : ℕ+) : ℕ)
          • (gd + ((ElemFrobCat.Hom.deg g : ℕ+) : ℕ) • c₁ + c₂)
        = gd + ((ElemFrobCat.Hom.deg g : ℕ+) : ℕ) • c₁
    rw [zmod2_odd_nsmul _ hodd]
    have h : c₂ + (gd + (((ElemFrobCat.Hom.deg g : ℕ+) : ℕ) • c₁) + c₂)
        = (c₂ + c₂) + (gd + (((ElemFrobCat.Hom.deg g : ℕ+) : ℕ) • c₁)) := by abel
    rw [h, zmod2_add_self, zero_add]
  · rw [ElemFrobCat.degFr_comp, ElemFrobCat.degFr_comp, z2_nfHom_deg, z2_nfHom_deg,
      z2Hom_deg]
    exact mul_comm _ _

theorem z2_naiveFrob_essSurj : (naiveFrob z2P z2F 3).EssSurj := by
  constructor
  intro B
  exact ⟨B, ⟨eqToIso (Subsingleton.elim _ _)⟩⟩

/-- ★★★**`Ψ_3` は `𝔽_{ℤ/2}` 上で圏同値**。 -/
theorem z2_naiveFrob_isEquivalence : (naiveFrob z2P z2F 3).IsEquivalence :=
  ⟨z2_naiveFrob_faithful, z2_naiveFrob_full, z2_naiveFrob_essSurj⟩

/-- ★`ℤ/2` 上の次数 2 の Frobenius 型射。 -/
def z2Frob2 : constObj (ZMod 2) ⟶ constObj (ZMod 2) := z2Hom _ _ 0 2

/-- ★`Div = 1` の pre-step。 -/
def z2Step1 : constObj (ZMod 2) ⟶ constObj (ZMod 2) := z2Hom _ _ 1 1

theorem z2Frob2_frobType : IsFrobeniusType z2P z2Frob2 := by
  refine (elemFrob_frobeniusType_iff Z2Phi cx2_totEpi (fun _ => isPreDivisorial_zmod2)
    z2Frob2).mpr ⟨?_, ?_⟩
  · show toChar (ElemFrobCat.Hom.div z2Frob2) = 0
    rw [show ElemFrobCat.Hom.div z2Frob2 = (0 : ZMod 2) from rfl]
    exact map_zero _
  · show IsIso (ElemFrobCat.Hom.base z2Frob2)
    show IsIso (eqToHom (Subsingleton.elim _ _))
    infer_instance

theorem z2Frob2_deg : z2P.degFr z2Frob2 = 2 := rfl

theorem z2Step1_preStep : IsPreStep z2P z2Step1 := by
  refine ⟨rfl, ?_⟩
  show IsIso (ElemFrobCat.Hom.base z2Step1)
  show IsIso (eqToHom (Subsingleton.elim _ _))
  infer_instance

/-- ★★★**`𝔽_{ℤ/2}` は perfect 型でない** —— `n = 2` で破れる
(`ℤ/2` では `2 • x = 0` なので `2 • a = 1` に解が無い)。 -/
theorem z2_not_perfectType : ¬ IsOfPerfectType z2P := by
  intro h
  obtain ⟨-, huniq⟩ := h (constObj (ZMod 2)) 2
  obtain ⟨ψ, ⟨-, hsq⟩, -⟩ := huniq _ _ _ _ z2Frob2 z2Frob2
    z2Frob2_frobType z2Frob2_deg z2Frob2_frobType z2Frob2_deg
    ⟨Iso.refl _⟩ ⟨Iso.refl _⟩ z2Step1 z2Step1_preStep
  obtain ⟨pd, hpd⟩ := z2_div_exists ψ
  have hd := congrArg ElemFrobCat.Hom.div hsq
  rw [ElemFrobCat.div_comp, ElemFrobCat.div_comp] at hd
  simp only [z2Frob2, z2Step1, constPhi_map, z2Hom_div, z2Hom_deg, hpd] at hd
  have hd2 : (1 : ZMod 2) + ((1 : ℕ+) : ℕ) • (0 : ZMod 2)
      = (0 : ZMod 2) + ((2 : ℕ+) : ℕ) • pd := hd
  rw [zmod2_even_nsmul _ (show ((2 : ℕ+) : ℕ) % 2 = 0 from rfl)] at hd2
  simp only [smul_zero, add_zero, zero_add] at hd2
  exact absurd hd2 (by decide)

/-- ★★★**`d = 3`(すなわち `d ≥ 2`)でも「`d` 固定」の読みは偽である**。

★`Ψ_3` は圏同値だが `𝔽_{ℤ/2}` は perfect 型でない。 -/
theorem prop_2_1_iii_fixed_degree_false_deg_three :
    (naiveFrob z2P z2F 3).IsEquivalence ∧ ¬ IsOfPerfectType z2P :=
  ⟨z2_naiveFrob_isEquivalence, z2_not_perfectType⟩

end ABC3.Check.FrdI
