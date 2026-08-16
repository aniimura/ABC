import ABC3.Found.FrdI.Def31
import ABC3.Found.FrdI.Prop22

/-!
# [FrdI] Remark 3.1.1 —— iso-subanchor は決して isotropic でない

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.57。

原文 (FrdI p.57):
> An iso-subanchor of the Frobenioid C is never isotropic. [In particular, if

## ★中身(測定)

**主張は 2 つ**:

| # | 内容 |
|---|---|
| 1 | `𝒞` の iso-subanchor は決して isotropic でない |
| 2 | ゆえに `𝒞` が isotropic 型なら quasi-isotropic 型 |

## ★原文の証明(3 段)

1. **`Proposition 1.10, (iv)`** —— isotropic な `A` からは、素数ごとに
   irreducible な射があり、次数が違えば非同型。★**素数は無限にあるので
   同型類は無限個**、すなわち `A` は anchor ではない
2. **`Definition 1.3, (vii), (b)`**(`isotropicClosed`)—— isotropic な対象からの
   射の終域も isotropic。★よって **subanchor も isotropic でない**
3. **mono-minimal 性** —— iso-subanchor `A` が isotropic だとすると、
   subanchor `B` からの mono-minimal categorical quotient `B ⟶ A` が
   isotropification(`Proposition 1.9, (v)`)を通って `B ⟶ B^istr ⟶ A` と分解し、
   `B ⟶ B^istr` は isotropic hull ゆえ mono。★**mono-minimal 性から
   `B ⟶ B^istr` が同型**、すなわち `B` が isotropic —— 2 に矛盾

★**3 段とも実装した。** 3 の急所は `G` の `B^istr` への作用の両立性で、
`Aut B → Aut B^istr` が単射準同型であること(`hullAutHom` / `hullAutHom_injective`)を
`Proposition 2.2, (iv)` から得る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★段 1 —— isotropic な対象は anchor でない

原文 (FrdI p.57):
> Indeed, by Proposition 1.10, (iv), an anchor is never isotropic. Thus, by Definition
-/

include P in
/-- ★★★**原文の段 1** —— isotropic な対象は **anchor ではない**。

★**理由は素数の無限性**: `Proposition 1.10, (iv)` が各素数 `p` について
次数 `p` の irreducible 射を与え、次数が違えば `_A𝒞` で非同型。
★anchor なら同型類が有限個なので、素数から有限集合への**単射**ができてしまう。 -/
theorem not_isAnchor_of_isotropic (F : FrobenioidCore P) {A : C} (hA : IsIsotropic P A) :
    ¬ IsAnchor C A := by
  rintro ⟨s, hs, hcov⟩
  obtain ⟨hex, hne⟩ := prop_1_10_iv_infinitely_many P F hA
  -- 各素数に irreducible 射を選ぶ
  choose B φ hirr hdeg using hex
  -- その同型類の代表を `s` の中に選ぶ
  choose X hXs hXiso using fun (p : ℕ+) (hp : Nat.Prime (p : ℕ)) =>
    hcov (B p hp) (φ p hp) (hirr p hp)
  haveI : Finite ↥s := hs.to_subtype
  haveI : Infinite {n : ℕ // n.Prime} := Nat.infinite_setOf_prime.to_subtype
  -- 素数 → `s` の写像は単射
  have hinj : Function.Injective
      (fun n : {n : ℕ // n.Prime} =>
        (⟨X ⟨n.1, n.2.pos⟩ n.2 , hXs ⟨n.1, n.2.pos⟩ n.2⟩ : ↥s)) := by
    rintro ⟨p, hp⟩ ⟨q, hq⟩ heq
    by_contra hpq
    have hpq' : (⟨p, hp.pos⟩ : ℕ+) ≠ ⟨q, hq.pos⟩ := by
      intro h
      exact hpq (Subtype.ext (congrArg PNat.val h))
    obtain ⟨e⟩ := hXiso ⟨p, hp.pos⟩ hp
    obtain ⟨e'⟩ := hXiso ⟨q, hq.pos⟩ hq
    have hXeq : X ⟨p, hp.pos⟩ hp = X ⟨q, hq.pos⟩ hq := congrArg Subtype.val heq
    -- `Under A` での同型を作る
    have ee : Under.mk (φ ⟨p, hp.pos⟩ hp) ≅ Under.mk (φ ⟨q, hq.pos⟩ hq) :=
      e ≪≫ (hXeq ▸ (Iso.refl _)) ≪≫ e'.symm
    -- 底の同型と四角形
    refine hne (B ⟨p, hp.pos⟩ hp) (B ⟨q, hq.pos⟩ hq)
      (φ ⟨p, hp.pos⟩ hp) (φ ⟨q, hq.pos⟩ hq) ?_ ((Under.forget A).mapIso ee) ?_
    · rw [hdeg, hdeg]
      exact hpq'
    · exact Under.w ee.hom
  exact absurd (Finite.of_injective _ hinj) (by
    simpa using (inferInstance : Infinite {n : ℕ // n.Prime}).not_finite)

include P in
/-- ★★★**原文の段 2** —— isotropic な対象は **subanchor でもない**。

★`Definition 1.3, (vii), (b)`(`isotropicClosed`)—— isotropic な対象からの
射の終域も isotropic —— を段 1 に繋ぐだけ。 -/
theorem not_isSubanchor_of_isotropic (F : FrobenioidCore P) {A : C} (hA : IsIsotropic P A) :
    ¬ IsSubanchor C A := by
  rintro ⟨B, ψ, hB⟩
  exact not_isAnchor_of_isotropic P F (F.isotropicClosed ψ hA) hB

/-! ## ★段 3 —— iso-subanchor も isotropic でない

★**急所は「isotropic hull が `Aut` の単射準同型を誘導する」こと。**
`Proposition 2.2, (iv)`(`hullOTriHom` / `hullOTriHom_injective`)が
`End` の準同型と単射性を与えるので、単元に制限すれば `Aut B →* Aut B^istr` になる。

★それを `IsMonoMinimalQuotient` の `G ≃* G'` として渡すと、
`B ⟶ B^istr` が同型になり、`B` が isotropic になってしまう。
-/

section HullAut

variable {A B : C} (ψ : A ⟶ B) (hψ : IsIsotropicHull P ψ)

/-- ★**isotropic hull が誘導する `Aut A →* Aut B`**。

★`Proposition 2.2, (iv)` の `hullOTriHom : End A →* End B` を単元へ制限したもの。 -/
noncomputable def hullAutHom : Aut A →* Aut B :=
  ((Aut.unitsEndEquivAut B).toMonoidHom.comp
      (Units.map (hullOTriHom P ψ hψ))).comp (Aut.unitsEndEquivAut A).symm.toMonoidHom

theorem hullAutHom_hom (γ : Aut A) :
    ((hullAutHom P ψ hψ γ).hom : B ⟶ B) = hullOTriHom P ψ hψ (γ.hom : End A) := rfl

/-- ★`hullAutHom` を特徴づける四角形。 -/
theorem hullAutHom_sq (γ : Aut A) :
    ((γ.hom : End A) : A ⟶ A) ≫ ψ = ψ ≫ ((hullAutHom P ψ hψ γ).hom : B ⟶ B) := by
  rw [hullAutHom_hom]
  exact hullOTriMap_sq P ψ hψ (γ.hom : End A)

include P in
/-- ★**`hullAutHom` は単射** —— `hullOTriHom_injective`(`Proposition 2.2, (iv)`)から。 -/
theorem hullAutHom_injective (F : FrobenioidCore P) :
    Function.Injective (hullAutHom P ψ hψ) := by
  intro x y h
  refine Aut.ext ?_
  exact hullOTriHom_injective P F ψ hψ
    (congrArg (fun z : Aut B => (z.hom : B ⟶ B)) h)

end HullAut

include P in
/-- ★★★**[FrdI] Remark 3.1.1 の第 1 主張** —— **iso-subanchor は決して isotropic でない**。

★段 2(subanchor は isotropic でない)に、mono-minimal 性で帰着させる。 -/
theorem not_isIsoSubanchor_of_isotropic (F : FrobenioidCore P) {A : C}
    (hA : IsIsotropic P A) : ¬ IsIsoSubanchor C A := by
  rintro ⟨B, hBsub, G, ψ, -, hmm⟩
  -- `B` の isotropic hull
  obtain ⟨Bi, ι, hι⟩ := F.isotropicHullExists B
  haveI : Mono ι := F.preStepMono ι hι.2.1
  -- `A` は isotropic なので `ψ` は hull を通る
  obtain ⟨ψ', hψ', -⟩ := hι.2.2.2 A hA ψ
  -- `G` を `Aut Bi` へ単射に写す
  have hinj := hullAutHom_injective P ι hι F
  have hcompat : ∀ γ : G,
      (((γ : Aut B).hom : End B) : B ⟶ B) ≫ ι
        = ι ≫ (((G.equivMapOfInjective (hullAutHom P ι hι) hinj γ : Aut Bi)).hom : Bi ⟶ Bi) := by
    intro γ
    rw [Subgroup.coe_equivMapOfInjective_apply]
    exact hullAutHom_sq P ι hι (γ : Aut B)
  -- mono-minimal 性で `ι` が同型
  haveI : IsIso ι :=
    hmm Bi ι ψ' inferInstance hψ' (G.map (hullAutHom P ι hι))
      (G.equivMapOfInjective (hullAutHom P ι hι) hinj) hcompat
  -- よって `B` は isotropic —— 段 2 に矛盾
  exact not_isSubanchor_of_isotropic P F (isIsotropic_of_iso P (asIso ι) hι.2.2.1) hBsub

include P in
/-- ★★★**[FrdI] Remark 3.1.1 の第 2 主張** —— `𝒞` が isotropic 型なら
**quasi-isotropic 型**。

★`isotropic 型`ならどの対象も isotropic なので `¬ IsIsotropic` は空虚に偽、
第 1 主張により `IsIsoSubanchor` も偽。両者が同値になる。 -/
theorem isOfQuasiIsotropicType_of_isOfIsotropicType (F : FrobenioidCore P)
    (h : IsOfIsotropicType P) : IsOfQuasiIsotropicType C P := by
  intro A
  constructor
  · intro hn
    exact absurd (h A) hn
  · intro hsub
    exact absurd (h A) (fun hA => not_isIsoSubanchor_of_isotropic P F hA hsub)

/-! ## ★★★出典の紐付け(`.src`) -/

/-- ★★★**[FrdI] Remark 3.1.1** —— 2 主張が実装された。

| # | 主張 | 実装 |
|---|---|---|
| 1 | iso-subanchor は決して isotropic でない | `not_isIsoSubanchor_of_isotropic`(＋ `not_isAnchor_of_isotropic` / `not_isSubanchor_of_isotropic`) |
| 2 | isotropic 型 ⟹ quasi-isotropic 型 | `isOfQuasiIsotropicType_of_isOfIsotropicType` |

★段 3 の急所は **isotropic hull が `Aut` の単射準同型を誘導すること**
(`hullAutHom` / `hullAutHom_injective`)—— `Proposition 2.2, (iv)` を使う。 -/
def not_isIsoSubanchor_of_isotropic.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 57, item := "Remark 3.1.1",
    sectionId := "frdi-remark-3-1-1" }

end ABC3.Found.FrdI
