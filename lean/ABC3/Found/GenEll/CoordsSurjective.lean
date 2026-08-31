/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.CommonGluedRatio
import ABC3.Found.GenEll.GlobalChartSurjective
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★段 E3c が閉じた —— 座標族を実際に作る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★これは何か —— 段 E3c の終点

> `f : X ⟶ Spec A` が有限型、`M` が有限アフィン被覆で自明化し、`X_s` がアフィンなら、
> ★**ある `n`・`N` と座標族 `sc : Fin (N+1) → sheafify(M^{⊗(n+1)})(⊤)`** があって
> `sc 0 = s^{⊗(n+1)}` かつ **`A⁰_{x_0} → Γ(X, X_{sc 0})` が全射**。

★★これが段 E3d（`§9-834`）の `IsClosedImmersion.of_surjective_of_isAffine` に
**そのまま渡せる形**である。

## ★★★★★★組み立て

| 段 | 出典 | 役割 |
|---|---|---|
| 有限生成 ⟹ 全射判定 | `§9-833` | 試験元 `T` を有限個に絞る |
| 試験元の族に単一の指数 | `§9-846` | `T` 全部に効く `n` と切断 `t_k` |
| 大域の比は像に入る | `§9-842` | `s_j/s_i` の形なら像 |
| 底環の像は自動 | `§9-843` | 定数は次数 0 の元 |
| ★添字の付け替え | 本ファイル `coordFamily` | `sd` と `t_k` を `Fin (T.card+1)` に並べる |
| ★★開の同一視の輸送 | 本ファイル `resEquivOfEq` | `X_{s^{⊗(n+1)}} = X_s` を環同型で運ぶ |

## ★測定の記録

★★★**循環しかけた**——試験元 `T` はチャート `X_{sc 0}` の座標環の生成元だが、
`sc` は `T` から作る。★しかし `X_{s^{⊗(n+1)}} = X_s` は **`n` に依らない**
（`nonVanishing_unit_secPow`、`§9-844`）ので、`T` を先に `Γ(X, X_s)` で取れば循環しない。
★★そこを環同型 `resEquivOfEq` で運ぶのが本ファイルの配管である。

## ★残っている段（明示）

★★★これで段 E3c・E3c-2 は閉じた。段 E の残りは **E1d**（チャートの射を貼って
`ψ : X ⟶ ℙᴺ_R` を作る——`§9-836` に機構はあるが「重なりで一致する」ことが未着手）である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}} {A : CommRingCat.{0}}

/-! ## ★添字の付け替え -/

/-- ★**分母 1 つと有限個の分子を `Fin (card+1)` に並べる**。 -/
noncomputable def coordFamily {κ : Type} (T : Finset κ) {B : Type} (sd : B) (t : κ → B) :
    Fin (T.card + 1) → B :=
  Fin.cases sd (fun j => t ((T.equivFin.symm j) : κ))

theorem coordFamily_zero {κ : Type} (T : Finset κ) {B : Type} (sd : B) (t : κ → B) :
    coordFamily T sd t 0 = sd := rfl

theorem coordFamily_succ {κ : Type} (T : Finset κ) {B : Type} (sd : B) (t : κ → B)
    (j : Fin T.card) : coordFamily T sd t j.succ = t ((T.equivFin.symm j) : κ) := rfl

/-- ★★**どの試験元も座標のどれかになっている**。 -/
theorem coordFamily_mem {κ : Type} (T : Finset κ) {B : Type} (sd : B) (t : κ → B)
    {k : κ} (hk : k ∈ T) : ∃ j : Fin (T.card + 1), coordFamily T sd t j = t k := by
  refine ⟨(T.equivFin ⟨k, hk⟩).succ, ?_⟩
  rw [coordFamily_succ, Equiv.symm_apply_apply]

/-! ## ★★開が等しいときの切断の輸送 -/

/-- ★★**開集合が等しいときの切断の環同型**。

★`§9-844` の `nonVanishing_unit_secPow`（`X_{s^{⊗n}} = X_s`）を渡るために要る
——2 つの開は等しいが**型は違う**からである。 -/
noncomputable def resEquivOfEq {A' B' : X.Opens} (h : A' = B') :
    (Γ(X, B') : Type) ≃+* (Γ(X, A') : Type) where
  toFun := X.presheaf.map (homOfLE (le_of_eq h)).op
  invFun := X.presheaf.map (homOfLE (le_of_eq h.symm)).op
  left_inv := fun x => by rw [res_trans, res_self]
  right_inv := fun x => by rw [res_trans, res_self]
  map_mul' := fun x y => by rw [map_mul]
  map_add' := fun x y => by rw [map_add]

/-- ★**同型で運んだ写像が全射なら元も全射である**。 -/
theorem surjective_of_comp_resEquiv {A' B' : X.Opens} (h : A' = B') {R : Type} [CommRing R]
    (ψ : R →+* (Γ(X, A') : Type))
    (hs : Function.Surjective ((resEquivOfEq h).symm.toRingHom.comp ψ)) :
    Function.Surjective ψ := by
  intro y
  obtain ⟨r, hr⟩ := hs ((resEquivOfEq h).symm y)
  refine ⟨r, ?_⟩
  have := congrArg (resEquivOfEq h) hr
  simpa using this

/-! ## ★★★★★★★★★★段 E3c の終点 -/

set_option maxHeartbeats 4000000 in
/-- ★★★★★★★★★★**段 E3c が閉じた** —— 座標族を実際に作る。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`f : X ⟶ Spec A` が有限型、`M` が有限アフィン被覆で自明化し、`X_s` がアフィンなら、
★**ある `n`・`N` と座標族 `sc : Fin (N+1) → sheafify(M^{⊗(n+1)})(⊤)`** があって

* `sc 0 = s^{⊗(n+1)}`（分母は与えられた切断の冪）
* ★★**`A⁰_{x_0} → Γ(X, X_{sc 0})` が全射**

★★★これが段 E3d（`§9-834`）に**そのまま渡せる形**である。

## ★測定の記録

★**循環しかけた**——試験元 `T` はチャート `X_{sc 0}` の座標環の生成元だが、
`sc` は `T` から作る。★★しかし `X_{s^{⊗(n+1)}} = X_s` は **`n` に依らない**ので、
`T` を先に `Γ(X, X_s)` で取れば循環しない。そこを `resEquivOfEq` で運ぶ。 -/
theorem exists_coords_surjective (f : X ⟶ Spec A) [LocallyOfFiniteType f]
    {ι : Type} [Fintype ι] (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (s : (M.obj (op ⊤) : Type)) (haff : IsAffineOpen (nonVanishing M s))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (hφ : ∀ r : (Γ(Spec A, (⊤ : (Spec A).Opens)) : Type),
      f.appLE ⊤ ⊤ (by simp) r ∈ Set.range φ) :
    ∃ (n N : ℕ) (sc : Fin (N + 1) →
        (((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val.obj
          (op (⊤ : X.Opens)) : Type)),
      sc 0 = ((sheafifyUnit X (presheafTensorPow M (n + 1))).app (op ⊤)).hom
              (secPow M s (n + 1))
      ∧ Function.Surjective (globalAwayHom
          ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
          (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow hM (n + 1)))
          φ sc 0) := by
  classical
  obtain ⟨T, hT⟩ := exists_finset_surjective_criterion f haff
  obtain ⟨n, t, hkt⟩ := exists_common_glued_globalRatio T M hM U hcov hU hUij e s id
  refine ⟨n, T.card, coordFamily T
    (((sheafifyUnit X (presheafTensorPow M (n + 1))).app (op ⊤)).hom (secPow M s (n + 1))) t,
    rfl, ?_⟩
  have hEq : nonVanishing ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
      (coordFamily T (((sheafifyUnit X (presheafTensorPow M (n + 1))).app (op ⊤)).hom
        (secPow M s (n + 1))) t 0)
      = nonVanishing M s := nonVanishing_unit_secPow M hM s n
  refine surjective_of_comp_resEquiv hEq _ ?_
  refine hT _ ?_ ?_
  · rintro _ ⟨r, rfl⟩
    obtain ⟨z, hz⟩ := range_appLE_subset_global f _
      (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow hM (n + 1)))
      φ _ 0 hφ ⟨r, rfl⟩
    refine ⟨z, ?_⟩
    show X.presheaf.map (homOfLE (le_of_eq hEq.symm)).op (globalAwayHom _ _ φ _ 0 z) = _
    rw [hz]
    exact appLE_res f ⊤ _ _ (by simp) (le_of_eq hEq.symm) r
  · rintro k hk
    have hkT : k ∈ T := hk
    obtain ⟨j, hj⟩ := coordFamily_mem T _ t hkT
    have hk' := hkt k hkT
    simp only [id_eq] at hk'
    obtain ⟨z, hz⟩ := mem_range_globalAwayHom _
      (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow hM (n + 1)))
      φ _ 0 j (X.presheaf.map (homOfLE (le_of_eq hEq)).op k)
      (by rw [hj]; exact hk')
    refine ⟨z, ?_⟩
    show X.presheaf.map (homOfLE (le_of_eq hEq.symm)).op (globalAwayHom _ _ φ _ 0 z) = _
    rw [hz, res_trans, res_self]

/-! ## ★出典の紐付け(`.src`) -/

def coordFamily.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(分母 1 つと有限個の分子を Fin (card+1) に並べる)",
    sectionId := "genell-prop-1-4" }

def resEquivOfEq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(開集合が等しいときの切断の環同型)",
    sectionId := "genell-prop-1-4" }

def exists_coords_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E3c が閉じた——座標族を実際に作る)",
    sectionId := "genell-prop-1-4" }

def exists_coords_surjective.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_finset_surjective_criterion(有限生成 ⟹ 全射判定、§9-833)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_finset_surjective_criterion") 2,
    .citation "[ABC3]" "exists_common_glued_globalRatio(試験元の族に単一の指数、§9-846)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_common_glued_globalRatio") 2,
    .citation "[ABC3]" "mem_range_globalAwayHom(§9-842) / range_appLE_subset_global(§9-843)"
      (.inProject "ABC3" "ABC3.Found.GenEll.mem_range_globalAwayHom") 2,
    .citation "[Stacks]" "Lemma 01PW(の消費側——チャートの座標環が生成される段)"
      (.absent "mathlib に ample は無い(2026-08-28 実測)。本定理は前層加群の言葉で独自に組んだ") 7,
    .implicitStep
      ("★**循環しかけた**——試験元 T はチャート X_{sc 0} の座標環の生成元だが sc は T から作る。" ++
       "★★しかし X_{s^{⊗(n+1)}} = X_s は **n に依らない**(§9-844)ので、" ++
       "T を先に Γ(X, X_s) で取れば循環しない。そこを resEquivOfEq で運ぶ") 5,
    .implicitStep
      ("★★★これで段 E3c・E3c-2 は閉じた。段 E の残りは **E1d**" ++
       "(チャートの射を貼って ψ : X ⟶ ℙᴺ_R を作る——§9-836 に機構はあるが" ++
       "「重なりで一致する」ことが未着手)である") 6 ]

end ABC3.Found.GenEll
