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

/-- ★★★**(ii) の後半を `𝒞` の言葉で** —— isotropic な対象は Frobenius-trivial。

★`prop_2_5_ii_frobTrivial` は `istrPre P F` の言葉で書かれている。
★★`Istr P` は**充満部分圏**なので `End` は一致し、`istrPre` の `Base`/`degFr` は
包含関手を通して定義されているから、**そのまま落ちる**。 -/
theorem isFrobeniusTrivial_of_isotropic (Fc : FrobenioidCore P)
    (hmt : ∀ A : C, IsMetricallyTrivial P A)
    (X : C) (hX : IsIsotropic P X) : IsFrobeniusTrivial P X := by
  obtain ⟨ζ, hdeg, hprop⟩ := prop_2_5_ii_frobTrivial P Fc hmt ⟨X, hX⟩
  refine ⟨⟨⟨fun n => ((ζ n).hom : X ⟶ X), ?_⟩, ?_⟩, ?_, ?_⟩
  · exact congrArg (fun z : End (⟨X, hX⟩ : Istr P) => (z.hom : X ⟶ X)) (map_one ζ)
  · intro a b
    exact congrArg (fun z : End (⟨X, hX⟩ : Istr P) => (z.hom : X ⟶ X)) (map_mul ζ a b)
  · intro n; exact hdeg n
  · -- ★co-angular だけは `istrPre` と `P` で意味が違う(分解を取る圏が違う)。
    -- ★★`X` は isotropic なので `Proposition 1.4, (i)` が `C` の側で供給する。
    intro n
    exact ⟨(hprop n).1, ⟨⟨prop_1_4_i P _ (fun Z g => Fc.isotropicClosed g hX),
      (hprop n).2.1.2⟩, (hprop n).2.2⟩⟩

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

/-! ## ★★isotropic hull に沿った `𝒪^▷` の移送 —— 非 isotropic への拡張の部品

原文 (FrdI p.49):
> [which applies even if A is not isotropic — cf. Definition 2.3, (a), (b)] to

★★**原文は `Proposition 2.5, (iii)` の証明で「分裂は `A` が isotropic でなくても
使える」と述べ、根拠に `Definition 2.3, (a), (b)` を挙げる。**
★我々はここを機械で追おうとして、**(a)(b) だけからは出ないこと**を見つけた
(下の記録を見よ)。ここではその途中まで——確実に出る 2 本——を用意する。
-/

section HullTransport

variable (F : FrobenioidCore P) {A B : C} {φ : A ⟶ B} (hφ : IsIsotropicHull P φ)

include P hφ in
/-- ★★**移送は `Div` を底に沿って戻す** ——
`Div α = Φ.map (Base φ) (Div (hullOTriHom α))`。

★`φ` は等長 pre-step(`Div φ = 0`, `degFr φ = 1`)で `α` は base-identity なので、
`hullOTriMap_sq` の四角形の両辺の `Div` を取るとこれだけが残る。 -/
theorem hullOTriHom_div (α : End A) (hα : α ∈ OTri P A) :
    P.Div ((α : A ⟶ A)) = Φ.map (P.Base φ) (P.Div ((hullOTriHom P φ hφ α : B ⟶ B))) := by
  have hsq := hullOTriMap_sq P φ hφ α
  have hl : P.Div ((α : A ⟶ A) ≫ φ) = P.Div ((α : A ⟶ A)) := by
    rw [P.Div_comp, show P.Div φ = 0 from hφ.1,
      show P.Base ((α : A ⟶ A)) = 𝟙 _ by
        have h : P.Base ((α : A ⟶ A)) = P.Base (𝟙 A) := hα.1
        rwa [P.Base_id] at h,
      MonoidOn.map_id, show P.degFr φ = 1 from hφ.2.1.1]
    simp
  have hr : P.Div (φ ≫ (hullOTriMap P φ hφ α : B ⟶ B))
      = Φ.map (P.Base φ) (P.Div ((hullOTriMap P φ hφ α : B ⟶ B))) := by
    rw [P.Div_comp, show P.Div φ = 0 from hφ.1]
    simp
  rw [← hl, hsq, hr]
  rfl

include P F hφ in
/-- ★★**`τ(A^istr)` は `𝒪^▷(A)` から一意に来る** ——
存在は `Definition 2.3, (b)`(`hullMem`)、一意性は
`Proposition 2.2, (iv)` の単射性。

★これが「非 isotropic な `A` に `τ(A)` を定める」ための道具である。 -/
theorem hullTau_existsUnique {τ : ∀ X : C, Submonoid (End X)}
    (hτ : IsCharacteristicSplitting P F τ) (t : End B) (ht : t ∈ τ B) :
    ∃! s : OTri P A, hullOTriHom P φ hφ ((s : End A)) = t := by
  obtain ⟨s, hs⟩ := hτ.hullMem hφ t ht
  refine ⟨s, hs, fun s' hs' => ?_⟩
  exact Subtype.ext (hullOTriHom_injective P F φ hφ (hs'.trans hs.symm))

end HullTransport

/-! ## ★★★訂正 —— 「非 isotropic でも使える」は導ける。私が見落としていた

★★一度「原文の主張は `Definition 2.3, (a), (b)` から導けない」と記録したが、
★★★**それは誤りだった**。見落としていたのは `Definition 1.3, (iii), (b)` である。

原文 (FrdI p.24):
> angular pre-step of C, then any morphism A

★**恒等射は co-angular な pre-step** なので、これを `A′ = A` に当てると
★★★**`A` の自己射はすべて co-angular**——**isotropy は要らない**。

★★**しかもこの補題は既に `Prop18.lean` にあった**(`isCoAngular_of_endo`)。
★自分たちで作った道具を使い損ねて「導けない」と記録しかけたのである。
★★**教訓: 「導けない」と書く前に、まず自分の在庫を検索する。**

★★これで `Proposition 2.2, (iii)`(`otri_div_eq_iff`、「`Div` が等しい ⟺
`𝒪^×` 倍だけ違う」)が**任意の `A` で成り立つ**(下の `otri_div_eq_iff'`)。
★詰まっていた「`u'` が像に入るか」は**迂回できる**——`u'` を持ち上げる代わりに、
`Div x = Div t` から直接 `x = t ≫ u` を得ればよい。

★**残る作業**(次段): `τ(A)` を非 isotropic な `A` で `τ(A^istr)` の引き戻しとして
確定させること。原文の `τ` は `(𝒞^istr)^lin` 上の部分関手なので、
非 isotropic な `A` での `τ(A)` は**その引き戻しが定義**である。
我々は `τ` を全対象で与えているため、この対応を条件として書く必要がある。
-/

include P in
/-- ★★★**`Proposition 2.2, (iii)` は任意の `A` で成り立つ** ——
`otri_div_eq_iff` の isotropy 仮定を外したもの。

★co-angular 性は `isCoAngular_of_endo` から来る。 -/
theorem otri_div_eq_iff' (F : FrobenioidCore P) {A : C} (x y : OTri P A) :
    P.Div ((x : End A) : A ⟶ A) = P.Div ((y : End A) : A ⟶ A)
      ↔ ∃ u : OTimes P A, ((x : End A) : A ⟶ A) = ((y : End A) : A ⟶ A)
          ≫ ((u : End A) : A ⟶ A) := by
  have hbx : P.Base ((x : End A) : A ⟶ A) = 𝟙 _ := by
    have h : P.Base ((x : End A) : A ⟶ A) = P.Base (𝟙 A) := x.2.1
    rwa [P.Base_id] at h
  have hby : P.Base ((y : End A) : A ⟶ A) = 𝟙 _ := by
    have h : P.Base ((y : End A) : A ⟶ A) = P.Base (𝟙 A) := y.2.1
    rwa [P.Base_id] at h
  have hsx : IsPreStep P ((x : End A) : A ⟶ A) :=
    ⟨x.2.2, by show IsIso (P.Base ((x : End A) : A ⟶ A)); rw [hbx]; infer_instance⟩
  have hsy : IsPreStep P ((y : End A) : A ⟶ A) :=
    ⟨y.2.2, by show IsIso (P.Base ((y : End A) : A ⟶ A)); rw [hby]; infer_instance⟩
  constructor
  · intro h
    obtain ⟨u, hu, he⟩ := F.faithfulUpToUnits ((x : End A) : A ⟶ A) ((y : End A) : A ⟶ A)
      (show P.Base _ = P.Base _ from by rw [hbx, hby]) h
      (isCoAngular_of_endo P F _) hsx (isCoAngular_of_endo P F _) hsy
    exact ⟨⟨u, hu⟩, he⟩
  · rintro ⟨u, he⟩
    haveI : IsIso ((u : End A) : A ⟶ A) := (CategoryTheory.isUnit_iff_isIso _).mp u.2.2
    rw [he, P.Div_comp, hby, MonoidOn.map_id,
      show P.Div ((u : End A) : A ⟶ A) = 0 from isIsometric_of_isIso P _,
      show P.degFr ((u : End A) : A ⟶ A) = 1 from u.2.1.2]
    simp

/-! ## ★★★★分裂は **`A` が isotropic でなくても**全単射

原文 (FrdI p.49):
> [which applies even if A is not isotropic — cf. Definition 2.3, (a), (b)] to

★★原文のこの一言を、我々の道具 4 つで閉じる:

| 道具 | 役割 |
|---|---|
| `hullPullback`(`Definition 2.3, (b)`) | `τ(A)` は `τ(A^istr)` の引き戻し |
| `hullTau_existsUnique` | `τ(A^istr)` の元は `𝒪^▷(A)` から**一意に**来る |
| `hullOTriHom_div` | `Div` が hull を通して底に沿って対応する |
| `otri_div_eq_iff'` | ★**`Div` が等しい ⟺ `𝒪^×` 倍だけ違う**(isotropy 不要) |

★★★**最後の 1 本が鍵**だった——`Definition 1.3, (iii), (b)` を恒等射に当てると
自己射はすべて co-angular なので、`Proposition 2.2, (iii)` に isotropy は要らない。
-/

theorem charSplitting_bijective_all (F : FrobenioidCore P)
    {τ : ∀ X : C, Submonoid (End X)} (hτ : IsCharacteristicSplitting P F τ) (A : C) :
    Function.Bijective
      (fun p : OTimes P A × τ A =>
        (⟨((p.1 : End A)) * ((p.2 : End A)),
          mul_mem (OTimes_le_OTri P A p.1.2) (hτ.le_otri A p.2.2)⟩ : OTri P A)) := by
  obtain ⟨B, φ, hφ⟩ := F.isotropicHullExists A
  have hunit : ∀ u : OTimes P A, P.Div (((u : End A)) : A ⟶ A) = 0 := by
    intro u
    haveI : IsIso (((u : End A)) : A ⟶ A) := (CategoryTheory.isUnit_iff_isIso _).mp u.2.2
    exact isIsometric_of_isIso P _
  have hdiv : ∀ (u : OTimes P A) (t : OTri P A),
      P.Div ((((u : End A)) * ((t : End A)) : End A) : A ⟶ A)
        = P.Div (((t : End A)) : A ⟶ A) := by
    intro u t
    show P.Div ((((t : End A)) : A ⟶ A) ≫ (((u : End A)) : A ⟶ A)) = _
    rw [P.Div_comp, hunit u,
      show P.Base (((t : End A)) : A ⟶ A) = 𝟙 _ from by
        have h : P.Base (((t : End A)) : A ⟶ A) = P.Base (𝟙 A) := t.2.1
        rwa [P.Base_id] at h,
      MonoidOn.map_id, show P.degFr (((u : End A)) : A ⟶ A) = 1 from u.2.1.2]
    simp
  -- ★hull を通した `Div` の対応
  have hpull : ∀ z : OTri P A, P.Div (((z : End A)) : A ⟶ A)
      = Φ.map (P.Base φ) (P.Div ((hullOTriHom P φ hφ (z : End A) : B ⟶ B))) :=
    fun z => hullOTriHom_div P hφ (z : End A) z.2
  constructor
  · rintro ⟨u₁, t₁⟩ ⟨u₂, t₂⟩ h
    have h' : (((u₁ : End A)) * ((t₁ : End A)) : End A)
        = (((u₂ : End A)) * ((t₂ : End A)) : End A) := congrArg Subtype.val h
    have hd : P.Div (((t₁ : End A)) : A ⟶ A) = P.Div (((t₂ : End A)) : A ⟶ A) := by
      rw [← hdiv u₁ ⟨(t₁ : End A), hτ.le_otri A t₁.2⟩,
        ← hdiv u₂ ⟨(t₂ : End A), hτ.le_otri A t₂.2⟩]
      exact congrArg (fun z : End A => P.Div (z : A ⟶ A)) h'
    -- ★hull へ送って `charBij` の一意性を使う
    have hdB : P.Div ((hullOTriHom P φ hφ (t₁ : End A) : B ⟶ B))
        = P.Div ((hullOTriHom P φ hφ (t₂ : End A) : B ⟶ B)) := by
      refine Φ.map_injective (P.Base φ) ?_
      rw [← hpull ⟨(t₁ : End A), hτ.le_otri A t₁.2⟩, ← hpull ⟨(t₂ : End A), hτ.le_otri A t₂.2⟩]
      exact hd
    have hm₁ : hullOTriHom P φ hφ (t₁ : End A) ∈ τ B :=
      (hτ.hullPullback hφ _ (hτ.le_otri A t₁.2)).mp t₁.2
    have hm₂ : hullOTriHom P φ hφ (t₂ : End A) ∈ τ B :=
      (hτ.hullPullback hφ _ (hτ.le_otri A t₂.2)).mp t₂.2
    obtain ⟨w, -, huniq⟩ := hτ.charBij B hφ.2.2.1
      ⟨hullOTriHom P φ hφ (t₁ : End A), hullOTriHom_mem P φ hφ _ (hτ.le_otri A t₁.2)⟩
    have hiota : hullOTriHom P φ hφ (t₁ : End A) = hullOTriHom P φ hφ (t₂ : End A) :=
      congrArg Subtype.val
        ((huniq ⟨_, hm₁⟩ rfl).trans (huniq ⟨_, hm₂⟩ hdB.symm).symm)
    have ht : t₁ = t₂ := Subtype.ext (hullOTriHom_injective P F φ hφ hiota)
    subst ht
    refine Prod.ext ?_ rfl
    refine Subtype.ext ?_
    haveI : Epi ((((t₁ : End A))) : A ⟶ A) := P.totEpiC _ _ _
    exact (cancel_epi ((((t₁ : End A))) : A ⟶ A)).mp h'
  · intro x
    -- ★`ι x` を `A^istr` で分裂させ、`τ` 成分を `A` へ引き戻す
    obtain ⟨t', ht', -⟩ := hτ.charBij B hφ.2.2.1
      ⟨hullOTriHom P φ hφ (x : End A), hullOTriHom_mem P φ hφ _ x.2⟩
    obtain ⟨s, hs⟩ := hτ.hullMem hφ (t' : End B) t'.2
    have hsτ : ((s : End A)) ∈ τ A :=
      (hτ.hullPullback hφ _ s.2).mpr (by rw [hs]; exact t'.2)
    have hds : P.Div (((s : End A)) : A ⟶ A) = P.Div (((x : End A)) : A ⟶ A) := by
      rw [hpull s, hpull x, hs, ht']
    obtain ⟨u, hu⟩ := (otri_div_eq_iff' P F x ⟨(s : End A), s.2⟩).mp hds.symm
    exact ⟨⟨u, ⟨(s : End A), hsτ⟩⟩, Subtype.ext hu.symm⟩

/-! ## ★★★★`τ` は**同型による共役**で保たれる —— isotropy 不要

★★`IsCharacteristicSplitting.map_mem` は `IsIsotropic` を要求する
(原文の `τ` が `(𝒞^istr)^lin` 上の部分関手だから)。
★★★**しかし「射が同型の場合」だけは isotropy 抜きで出る。**

★筋: 同型 `j : Y ⟶ Y'` は hull を hull へ写す(`j ≫ h'` は `Y` の hull)。
hull の普遍性で `θ : H ⟶ H'` を取ると、`h` が mono なので
**共役が hull を通して対応する**:

  `θ ≫ (h' を通した t) = (h を通した s) ≫ θ`

★あとは **`H`・`H'` は isotropic** なので `map_mem` が使え、
`hullPullback` で `Y` 側へ引き戻す。

★★**これが `isPsiValue_unique`(well-defined 性)の要**である——
そこに現れる射は `arbFactorUniq` の出す**同型**だけだから。
-/

theorem tau_conj_mem (F : FrobenioidCore P) {τ : ∀ X : C, Submonoid (End X)}
    (hτ : IsCharacteristicSplitting P F τ) {Y Y' : C} (j : Y ⟶ Y') [IsIso j]
    (t : OTri P Y') (ht : ((t : End Y')) ∈ τ Y') :
    ((otriPull P F j (isCoAngular_of_isIso P j) (isLinear_of_isIso P j) t : End Y))
      ∈ τ Y := by
  obtain ⟨H, h, hh⟩ := F.isotropicHullExists Y
  obtain ⟨H', h', hh'⟩ := F.isotropicHullExists Y'
  -- ★手 1: hull の普遍性で `θ : H ⟶ H'`
  obtain ⟨θ, hθ, -⟩ := hh.2.2.2 H' hh'.2.2.1 (j ≫ h')
  have hθlin : IsLinear P θ := by
    have hd := congrArg P.degFr hθ
    rw [P.degFr_comp, P.degFr_comp, show P.degFr j = 1 from isLinear_of_isIso P j,
      show P.degFr h' = 1 from hh'.2.1.1, show P.degFr h = 1 from hh.2.1.1] at hd
    show P.degFr θ = 1
    simpa using hd.symm
  -- ★手 2: 共役が hull を通して対応する
  set s := otriPull P F j (isCoAngular_of_isIso P j) (isLinear_of_isIso P j) t with hs
  have hspec : j ≫ ((t : End Y') : Y' ⟶ Y') = ((s : End Y) : Y ⟶ Y) ≫ j :=
    otriPull_spec P F j (isCoAngular_of_isIso P j) (isLinear_of_isIso P j) t
  haveI : Epi h := P.totEpiC _ _ _
  have hsq : θ ≫ (hullOTriMap P h' hh' ((t : End Y')) : H' ⟶ H')
      = (hullOTriMap P h hh ((s : End Y)) : H ⟶ H) ≫ θ := by
    refine (cancel_epi h).mp ?_
    calc h ≫ θ ≫ (hullOTriMap P h' hh' ((t : End Y')) : H' ⟶ H')
        = (h ≫ θ) ≫ (hullOTriMap P h' hh' ((t : End Y')) : H' ⟶ H') := by simp
      _ = (j ≫ h') ≫ (hullOTriMap P h' hh' ((t : End Y')) : H' ⟶ H') := by rw [hθ]
      _ = j ≫ (h' ≫ (hullOTriMap P h' hh' ((t : End Y')) : H' ⟶ H')) := by simp
      _ = j ≫ (((t : End Y') : Y' ⟶ Y') ≫ h') := by rw [← hullOTriMap_sq P h' hh' _]
      _ = (j ≫ ((t : End Y') : Y' ⟶ Y')) ≫ h' := by simp
      _ = (((s : End Y) : Y ⟶ Y) ≫ j) ≫ h' := by rw [hspec]
      _ = ((s : End Y) : Y ⟶ Y) ≫ (j ≫ h') := by simp
      _ = ((s : End Y) : Y ⟶ Y) ≫ (h ≫ θ) := by rw [hθ]
      _ = (((s : End Y) : Y ⟶ Y) ≫ h) ≫ θ := by simp
      _ = (h ≫ (hullOTriMap P h hh ((s : End Y)) : H ⟶ H)) ≫ θ := by
            rw [hullOTriMap_sq P h hh _]
      _ = h ≫ (hullOTriMap P h hh ((s : End Y)) : H ⟶ H) ≫ θ := by simp
  -- ★手 3: `H` 側で `map_mem`、`Y` 側へ `hullPullback`
  have htH' : (hullOTriHom P h' hh' ((t : End Y'))) ∈ τ H' :=
    (hτ.hullPullback hh' _ t.2).mp ht
  have hmem := hτ.map_mem hh.2.2.1 hθlin
    (⟨hullOTriHom P h' hh' ((t : End Y')), hullOTriHom_mem P h' hh' _ t.2⟩ : OTri P H') htH'
  have heq : (⟨hullOTriHom P h hh ((s : End Y)),
      hullOTriHom_mem P h hh _ s.2⟩ : OTri P H)
      = otriLin P F hh.2.2.1 hθlin
        ⟨hullOTriHom P h' hh' ((t : End Y')), hullOTriHom_mem P h' hh' _ t.2⟩ :=
    otriLin_uniq P F hh.2.2.1 hθlin _ _ hsq
  refine (hτ.hullPullback hh _ s.2).mpr ?_
  rw [show hullOTriHom P h hh ((s : End Y))
    = ((otriLin P F hh.2.2.1 hθlin
        ⟨hullOTriHom P h' hh' ((t : End Y')), hullOTriHom_mem P h' hh' _ t.2⟩ :
          OTri P H) : End H) from congrArg Subtype.val heq]
  exact hmem

/-! ## ★★同型に沿った共役準同型 —— `isPsiValue_unique` が要る形

★`isPsiValue_unique`(`Ψ` の well-defined 性)に現れる射は
`arbFactorUniq` の出す**同型だけ**である。★したがってそこで使う道具を
**同型専用**に用意すれば、isotropy を仮定せずに済む。
-/

section Conj

variable (F : FrobenioidCore P) {Y Y' : C} (j : Y ⟶ Y') [IsIso j]

/-- ★★**同型に沿った共役** `𝒪^▷(Y') →* 𝒪^▷(Y)` —— isotropy 不要。

★同型は co-angular かつ linear なので `Proposition 1.11, (iv)` がそのまま使える。 -/
noncomputable def conjHom : OTri P Y' →* OTri P Y :=
  otriPullHom P F j (isCoAngular_of_isIso P j) (isLinear_of_isIso P j)

include P in
/-- ★共役の定義式 —— `j ≫ t = conj(t) ≫ j`。 -/
theorem conjHom_spec (t : OTri P Y') :
    j ≫ (((t : End Y')) : Y' ⟶ Y') = (((conjHom P F j t : OTri P Y) : End Y) : Y ⟶ Y) ≫ j :=
  otriPull_spec P F j (isCoAngular_of_isIso P j) (isLinear_of_isIso P j) t

include P in
/-- ★★**共役は単元を単元へ写す** —— モノイド準同型なので、
`𝒪^▷` の中での逆元がそのまま移る。 -/
theorem conjHom_otimes_mem (u : OTri P Y') (hu : ((u : End Y')) ∈ OTimes P Y') :
    ((conjHom P F j u : OTri P Y) : End Y) ∈ OTimes P Y := by
  haveI : IsIso ((((u : End Y'))) : Y' ⟶ Y') := (CategoryTheory.isUnit_iff_isIso _).mp hu.2
  have hbu : P.Base ((((u : End Y'))) : Y' ⟶ Y') = 𝟙 _ := by
    have h : P.Base ((((u : End Y'))) : Y' ⟶ Y') = P.Base (𝟙 Y') := u.2.1
    rwa [P.Base_id] at h
  have hinv : (inv ((((u : End Y'))) : Y' ⟶ Y')) ∈ OTri P Y' := by
    refine ⟨?_, degFr_of_isIso P _⟩
    show P.Base (inv ((((u : End Y'))) : Y' ⟶ Y')) = P.Base (𝟙 Y')
    have h : P.Base (inv ((((u : End Y'))) : Y' ⟶ Y'))
        ≫ P.Base ((((u : End Y'))) : Y' ⟶ Y') = P.Base (𝟙 Y') := by
      rw [← P.Base_comp, IsIso.inv_hom_id]
    rwa [hbu, Category.comp_id] at h
  have hmul : u * (⟨inv ((((u : End Y'))) : Y' ⟶ Y'), hinv⟩ : OTri P Y') = 1 :=
    Subtype.ext (by show inv ((((u : End Y'))) : Y' ⟶ Y') ≫ _ = _; simp)
  have hmul' : (⟨inv ((((u : End Y'))) : Y' ⟶ Y'), hinv⟩ : OTri P Y') * u = 1 :=
    Subtype.ext (by show ((((u : End Y'))) : Y' ⟶ Y') ≫ _ = _; simp)
  refine ⟨(conjHom P F j u).2, (CategoryTheory.isUnit_iff_isIso _).mpr ?_⟩
  refine ⟨((conjHom P F j ⟨inv ((((u : End Y'))) : Y' ⟶ Y'), hinv⟩ : OTri P Y) : End Y),
    ?_, ?_⟩
  · have h := congrArg (fun z : OTri P Y => ((z : End Y)))
      ((conjHom P F j).map_mul ⟨inv ((((u : End Y'))) : Y' ⟶ Y'), hinv⟩ u).symm
    simp only [hmul', map_one] at h
    exact h
  · have h := congrArg (fun z : OTri P Y => ((z : End Y)))
      ((conjHom P F j).map_mul u ⟨inv ((((u : End Y'))) : Y' ⟶ Y'), hinv⟩).symm
    simp only [hmul, map_one] at h
    exact h

include P in
/-- ★★★**共役は `τ` を保つ** —— `tau_conj_mem` の言い換え。 -/
theorem conjHom_tau_mem {τ : ∀ X : C, Submonoid (End X)}
    (hτ : IsCharacteristicSplitting P F τ) (t : OTri P Y')
    (ht : ((t : End Y')) ∈ τ Y') : ((conjHom P F j t : OTri P Y) : End Y) ∈ τ Y :=
  tau_conj_mem P F hτ j t ht

end Conj

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
