import ABC3.Found.FrdI.Def23

/-!
# [FrdI] Proposition 2.5 —— The Unit-linear Frobenius Functor

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.48–p.50。

原文 (FrdI p.48):
> (The Unit-linear Frobenius Functor) Let

## ★この命題の規模(測定)

**3 条、主張は 4**:

| 条 | # | 内容 |
|---|---|---|
| (i) | 1 | `𝒪^▷(A)^char → Φ(A)` が**全単射**(`Proposition 2.2, (iii)` の包含が実は全単射) |
| (ii) | 2 | `𝒞^istr` が **base-trivial 型** |
| (ii) | 3 | `𝒞^istr` のどの対象も **Frobenius-trivial** |
| (iii) | 4 | **unit-linear Frobenius 函手** `Ψ : 𝒞 ≃ 𝒞(d)`((a) 対象と等長射の上で恒等、(b) `d` の Frobenius 函手と 1-compatible) |

★**(iii) は `Definition 2.4` の `Λ` / `d ∈ Λ>0` / `𝒞(d)` を要する**ので、
ここでは (i)(ii) を実装する。★**(iii) が入るまで `.src` は付けない。**

## ★仮定について

原文は `𝒞` が **Frobenius-normalized・metrically trivial・Aut-ample 型**であるとし、
`τ` を characteristic splitting、`Λ` を `Φ` を supports する monoid type、`d ∈ Λ>0` とする。
★**(i) が使うのは metrically trivial と Aut-ample だけ**(原文の証明どおり)、
★**(ii) が使うのは metrically trivial だけ**である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★(i) —— `𝒪^▷(A)^char → Φ(A)` は全単射

原文 (FrdI p.48):
> sition 2.2, (iii), is, in fact, a bijection.

原文 (FrdI p.49):
> metrically trivial and Aut-ample type. Next, we consider assertion (ii). Since C

★**原文の証明どおりの 3 手**:
1. `Definition 1.3, (iii), (d)` の圏同値で、`c ∈ Φ(A)` を
   co-angular pre-step `ψ : A ⟶ X`(`Div ψ = c`)として実現する(`coaPre_realize`)
2. ★**metrically trivial 型**なので `X ≅ A`、すなわち `ψ` は `A` の自己射に化ける
3. ★**Aut-ample 型**なので、ずれた底の自己同型を `𝒞` の自己同型で打ち消せる

★**単射性の側は `Proposition 2.2, (iii)`(`otri_div_eq_iff`)で既に取れている。**
-/

/-- ★★★**[FrdI] Proposition 2.5, (i)** —— `Div : 𝒪^▷(A) → Φ(A)` は**全射**。

`Proposition 2.2, (iii)` の単射性(`otri_div_eq_iff`)と合わせて、
`𝒪^▷(A)^char → Φ(A)` が**全単射**であることを与える。 -/
theorem prop_2_5_i_surjective (G : Frobenioid P) {A : C}
    (hmt : IsMetricallyTrivial P A) (haa : IsAutAmple P A)
    (c : Φ.val (P.toElem.obj A).base) :
    ∃ x : OTri P A, P.Div ((x : End A) : A ⟶ A) = c := by
  -- 手 1: `c` を co-angular pre-step として実現する
  obtain ⟨X, ψ, hψc, hψs, hψd⟩ := coaPre_realize P G A c
  -- 手 2: metrically trivial 型で `X ≅ A`
  obtain ⟨e⟩ := hmt X ψ hψc hψs
  -- 手 3: Aut-ample 型でずれた底を打ち消す
  have hgi : IsIso (P.Base (ψ ≫ e.hom)) := by
    haveI : IsIso (P.Base ψ) := hψs.2
    haveI : IsIso (P.Base e.hom) := isBaseIsomorphism_of_isIso P e.hom
    rw [P.Base_comp]
    infer_instance
  obtain ⟨φ, hφi, hφb⟩ := haa (P.Base (ψ ≫ e.hom)) hgi
  haveI := hφi
  have hbinv : P.Base (ψ ≫ e.hom) ≫ P.Base (inv φ) = 𝟙 _ := by
    rw [← hφb, ← P.Base_comp, IsIso.hom_inv_id, P.Base_id]
  refine ⟨⟨ψ ≫ e.hom ≫ inv φ, ?_, ?_⟩, ?_⟩
  · -- base-identity
    show P.Base (ψ ≫ e.hom ≫ inv φ) = P.Base (𝟙 A)
    rw [P.Base_id, ← Category.assoc, P.Base_comp, hbinv]
  · -- linear
    show P.degFr (ψ ≫ e.hom ≫ inv φ) = 1
    rw [P.degFr_comp, P.degFr_comp, show P.degFr ψ = 1 from hψs.1,
      degFr_of_isIso P e.hom, degFr_of_isIso P (inv φ)]
    simp
  · -- `Div = c`
    show P.Div (ψ ≫ e.hom ≫ inv φ) = c
    haveI : IsIso (e.hom ≫ inv φ) := inferInstance
    rw [P.Div_comp, show P.Div (e.hom ≫ inv φ) = 0 from isIsometric_of_isIso P _,
      degFr_of_isIso P (e.hom ≫ inv φ), map_zero]
    simpa using hψd

/-! ## ★(ii) —— `𝒞^istr` は base-trivial 型で、どの対象も Frobenius-trivial

原文 (FrdI p.48):
> (ii) Cistr is of base-trivial type. Moreover, every object of Cistr is Frobenius-

原文 (FrdI p.49):
> follows from the existence of Frobenius-trivial objects [cf. Definition 1.3, (i), (a)]

★**原文の証明どおり**:
1. `Definition 1.3, (i), (b)`(`preStepSpan`)で、底が同型な `𝒞^istr` の 2 対象は
   pre-step の span で結ばれる。★**span の頂点も `𝒞^istr` の対象**なので
   `Proposition 1.4, (i)` により両辺の pre-step は自動で co-angular
2. ★**metrically trivial 型**なので、その span の両端は頂点と同型 —— よって互いに同型
3. `Definition 1.3, (i), (a)` を **`𝒞^istr` 自身**に当てる
   (`Proposition 1.9, (v)`、`istr_frobenioidCore`)と Frobenius-trivial 対象が取れ、
   2 の base-trivial 性で**すべての対象へ移る**
-/

variable (F : FrobenioidCore P)

/-- ★★★**[FrdI] Proposition 2.5, (ii) の前半** —— `𝒞^istr` は **base-trivial 型**。 -/
theorem prop_2_5_ii_baseTrivial (hmt : ∀ A : C, IsMetricallyTrivial P A) (A : Istr P) :
    IsBaseTrivial (istrPre P F) A := by
  intro Dd hbi
  obtain ⟨α⟩ := hbi
  obtain ⟨X, σ, τ, hσ, hτ, -⟩ :=
    (istr_frobenioidCore P F).preStepSpan A Dd α.hom inferInstance
  obtain ⟨eA⟩ := hmt X.obj A.obj σ.hom
    (isCoAngular_of_isotropic_dom P F X.property σ.hom) hσ
  obtain ⟨eD⟩ := hmt X.obj Dd.obj τ.hom
    (isCoAngular_of_isotropic_dom P F X.property τ.hom) hτ
  exact ⟨ObjectProperty.isoMk _ (eD ≪≫ eA.symm)⟩

/-- ★★★**[FrdI] Proposition 2.5, (ii) の後半** —— `𝒞^istr` のどの対象も
**Frobenius-trivial**。

★`Definition 1.3, (i), (a)` を `𝒞^istr` 自身に当て、前半の base-trivial 性で移す。 -/
theorem prop_2_5_ii_frobTrivial (hmt : ∀ A : C, IsMetricallyTrivial P A) (A : Istr P) :
    IsFrobeniusTrivial (istrPre P F) A := by
  obtain ⟨B, hB, ⟨e⟩⟩ := (istr_frobenioidCore P F).baseSurj
    ((istrPre P F).toElem.obj A).base
  obtain ⟨θ⟩ := prop_2_5_ii_baseTrivial P F hmt A B ⟨e.symm⟩
  exact isFrobeniusTrivial_of_iso (istrPre P F) (istr_frobenioidCore P F) θ hB

/-! ## ★★divisorial なモノイドでは `d` 倍が単射

★★**`Ψ` が忠実・充満であることの土台**である。原文は `Proposition 2.5, (iii)` の
最後で「`Ψ`・`𝒞(d)` の定義から明らか」と書くが、そこには
**「`t^d` から `t` が一意に戻る」**という事実が要る。

★★★**証明は `Gp M` を経由する**:
`d • x = d • y` なら `z := x − y ∈ Gp M` が `d • z = 0` を満たす。
★`0` は `M` の像なので **saturated** から `z` も `−z` も `M` の像であり、
その和が `0` なので **integral** ＋ **sharp** で `z = 0`。
-/

/-- ★`toGp` は零を保つ。 -/
theorem toGp_zero (M : Type w) [AddCommMonoid M] : toGp M 0 = 0 :=
  AddLocalization.mk_zero

theorem toGp_add {M : Type w} [AddCommMonoid M] (a b : M) :
    toGp M (a + b) = toGp M a + toGp M b := by
  show AddLocalization.mk (a + b) (0 : (⊤ : AddSubmonoid M)) = _
  rw [show (toGp M a + toGp M b) = AddLocalization.mk (a + b) (0 + 0 : (⊤ : AddSubmonoid M))
    from AddLocalization.mk_add a b 0 0, zero_add]

theorem toGp_nsmul {M : Type w} [AddCommMonoid M] (n : ℕ) (a : M) :
    toGp M (n • a) = n • toGp M a := by
  induction n with
  | zero => simpa using toGp_zero M
  | succ k ih => rw [succ_nsmul, toGp_add, ih, succ_nsmul]

/-- ★★★**divisorial なモノイドでは `n` 倍(`n ≥ 1`)は単射**。

★integral(`Gp` への単射)・saturated・sharp の 3 つがちょうど使われる。 -/
theorem nsmul_injective_of_divisorial {M : Type w} [AddCommMonoid M]
    (hint : IsIntegralMonoid M) (hsat : IsSaturatedMonoid M) (hsh : IsSharp M)
    {n : ℕ} (hn : 0 < n) {x y : M} (h : n • x = n • y) : x = y := by
  have hmem : ∀ w : Gp M, n • w = 0 → w ∈ Set.range (toGp M) := by
    intro w hw
    refine hsat w n hn ?_
    rw [hw]
    exact ⟨0, toGp_zero M⟩
  have hnz : n • (toGp M x - toGp M y) = 0 := by
    rw [smul_sub, ← toGp_nsmul, ← toGp_nsmul, h, sub_self]
  obtain ⟨a, ha⟩ := hmem _ hnz
  obtain ⟨b, hb⟩ := hmem (-(toGp M x - toGp M y)) (by rw [smul_neg, hnz, neg_zero])
  have hab : a + b = 0 := by
    refine hint ?_
    rw [toGp_add, ha, hb, add_neg_cancel, toGp_zero]
  have ha0 : a = 0 := hsh a ⟨⟨a, b, hab, by rw [add_comm]; exact hab⟩, rfl⟩
  refine hint ?_
  rw [← sub_eq_zero, ← ha, ha0, toGp_zero]

/-! ## ★★★`𝒪^▷(A)^char → Φ(A)` の全単射 —— (i) を閉じる

原文 (FrdI p.48):
> (i) The natural inclusion O

★**特性モノイド `𝒪^▷(A)^char` を「`Div` の核による合同」として作る。**
`Proposition 2.2, (iii)`(`otri_div_eq_iff`)が
「`Div` が等しい ⟺ `𝒪^×` 倍だけ違う」を与えるので、
★**この商はちょうど原文の `𝒪^▷(A)^char`(＝ `𝒪^▷(A)/𝒪^×(A)`)である**
(`otriChar_rel_iff` で確認する)。

★この作り方だと**単射性は構成から自動**になり、
残るのは `prop_2_5_i_surjective` の全射性だけになる。
-/

section Char

variable {A : C}

include P in
/-- ★**`𝒪^▷(A)` 上では `Div` は加法的** —— base-identity かつ linear なので
合成則 `Div (ψ ≫ φ) = Φ.map (Base ψ) (Div φ) + (degFr φ) • Div ψ` が
ただの和に潰れる。 -/
theorem otri_div_mul (x y : OTri P A) :
    P.Div ((((x * y : OTri P A) : End A)) : A ⟶ A)
      = P.Div ((x : End A) : A ⟶ A) + P.Div ((y : End A) : A ⟶ A) := by
  have hby : P.Base ((y : End A) : A ⟶ A) = 𝟙 _ := by
    have h : P.Base ((y : End A) : A ⟶ A) = P.Base (𝟙 A) := y.2.1
    rwa [P.Base_id] at h
  show P.Div (((y : End A) : A ⟶ A) ≫ ((x : End A) : A ⟶ A)) = _
  rw [P.Div_comp, hby, MonoidOn.map_id, show P.degFr ((x : End A) : A ⟶ A) = 1 from x.2.2]
  simp

/-- ★**`Div` をモノイド準同型として見る** —— 終域は `Φ(A)` を乗法的に見たもの。 -/
def otriDivHom (A : C) : OTri P A →* Multiplicative (Φ.val (P.toElem.obj A).base) where
  toFun x := Multiplicative.ofAdd (P.Div ((x : End A) : A ⟶ A))
  map_one' := by
    show Multiplicative.ofAdd (P.Div (𝟙 A)) = 1
    rw [P.Div_id]
    rfl
  map_mul' x y := by
    show Multiplicative.ofAdd (P.Div (((x * y : OTri P A) : End A) : A ⟶ A)) = _
    rw [otri_div_mul P x y]
    rfl

/-- ★★**特性モノイド `𝒪^▷(A)^char`** —— `Div` の核による商。 -/
def OTriChar (A : C) : Type v2 := (Con.ker (otriDivHom P A)).Quotient

noncomputable instance : Monoid (OTriChar P A) :=
  inferInstanceAs (Monoid (Con.ker (otriDivHom P A)).Quotient)

/-- ★`𝒪^▷(A) ↠ 𝒪^▷(A)^char`。 -/
def toOTriChar : OTri P A →* OTriChar P A := (Con.ker (otriDivHom P A)).mk'

/-- ★★`𝒪^▷(A)^char → Φ(A)` —— 原文の「natural inclusion」。 -/
def otriCharDiv : OTriChar P A →* Multiplicative (Φ.val (P.toElem.obj A).base) :=
  Con.kerLift (otriDivHom P A)

/-- ★**構成から単射**。 -/
theorem otriCharDiv_injective : Function.Injective (otriCharDiv P (A := A)) :=
  Con.kerLift_injective _

include P in
/-- ★★**この商がちょうど原文の `𝒪^▷(A)/𝒪^×(A)` である**ことの確認。

★`Proposition 2.2, (iii)`(`otri_div_eq_iff`)そのもの。 -/
theorem otriChar_rel_iff (F : FrobenioidCore P) (hiso : IsIsotropic P A) (x y : OTri P A) :
    toOTriChar P x = toOTriChar P y
      ↔ ∃ u : OTimes P A, ((x : End A) : A ⟶ A)
          = ((y : End A) : A ⟶ A) ≫ ((u : End A) : A ⟶ A) := by
  rw [← otri_div_eq_iff P F hiso x y]
  constructor
  · intro h
    exact (Con.eq _).mp h
  · intro h
    exact (Con.eq _).mpr h

include P in
/-- ★★★**[FrdI] Proposition 2.5, (i)** —— `𝒪^▷(A)^char → Φ(A)` は**全単射**。

★単射性は構成から、全射性は `prop_2_5_i_surjective` から。 -/
theorem prop_2_5_i_bijective (G : Frobenioid P)
    (hmt : IsMetricallyTrivial P A) (haa : IsAutAmple P A) :
    Function.Bijective (otriCharDiv P (A := A)) := by
  refine ⟨otriCharDiv_injective P, ?_⟩
  intro c
  obtain ⟨x, hx⟩ := prop_2_5_i_surjective P G hmt haa (Multiplicative.toAdd c)
  refine ⟨toOTriChar P x, ?_⟩
  show Multiplicative.ofAdd (P.Div ((x : End A) : A ⟶ A)) = c
  rw [hx]
  rfl

include P in
/-- ★**同型として述べたもの**。 -/
noncomputable def prop_2_5_i_equiv (G : Frobenioid P)
    (hmt : IsMetricallyTrivial P A) (haa : IsAutAmple P A) :
    OTriChar P A ≃* Multiplicative (Φ.val (P.toElem.obj A).base) :=
  MulEquiv.ofBijective _ (prop_2_5_i_bijective P G hmt haa)

end Char

end ABC3.Found.FrdI
