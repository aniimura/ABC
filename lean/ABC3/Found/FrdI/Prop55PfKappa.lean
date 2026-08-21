/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32Perfect

/-!
# [FrdI] `𝒞^pf` の根の標準射 —— 選択された同型の整合性を消す

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> tween the respective sets of morphisms between the images of two given objects of C

## ★★何のためのファイルか

`Proposition 5.5, (ii)` の `birat` の側(`Prop55BiratPf.lean`)では、
`α : A ⟶ A₁` が次数 `k` の Frobenius 型のとき `⟨A,1⟩ ≅ ⟨A₁,k⟩` を使う。
★★しかしその同型を `frobDegUniq` の `.choose` で取ると、
**`W` を大きくしたときの整合性が言えない**(自己同型だけずれ得る)。

★★★**逃げ道**: `𝒞^pf` の Frobenius 型射は **mono** である
(`pfRoot_frobTypeMono`、在庫)。★したがって同型を**方程式で特徴づければ一意**になる:

    `e ≫ κ_{A₁,k} = [α]`      (`κ_{A₁,k} : ⟨A₁,k⟩ ⟶ ⟨A₁,1⟩` は標準の次数 `k` Frobenius)

★存在は `pfRoot_frob_div`(在庫、**共通の終域への Frobenius 型射による割り算**)が与え、
同型であることは `κ` と `[α]` がどちらも mono であることから出る。
★★整合性 `rootLift (α ≫ a) = rootLift α ≫ rootStep a` は
**両辺に `κ` を合成して比べるだけ**で出る。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `pfKappa` | ★標準の次数 `k` の Frobenius 型射 `⟨A,k⟩ ⟶ ⟨A,1⟩` |
| `pfKappa_degFr` / `pfKappa_frobType` / `pfKappa_mono` | その 3 性質 |
| `rootLift` / `rootLift_spec` / `rootLift_isIso` | ★★★`⟨A,1⟩ ≅ ⟨A₁,k⟩`(方程式つき) |
| `rootStep` / `rootStep_spec` / `rootStep_isIso` | ★★★`⟨A₁,k⟩ ≅ ⟨A₁′,k·d⟩`(方程式つき) |
| `rootLift_comp` | ★★★★★★**整合性**(これが選択の問題を消す) |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-! ## ★1. 標準の次数 `k` の Frobenius 型射 -/

/-- ★`rtObj A 1 ⟶ A`(標準の同型 `rtExt A 1` の逆)。 -/
noncomputable def rtOneInv (A : C) : rtObj P F A 1 ⟶ A :=
  @inv _ _ _ _ (rtExt P F A 1) (isIso_rtExt_one P F A)

instance rtOneInv_isIso (A : C) : IsIso (rtOneInv (P := P) (F := F) A) :=
  @IsIso.inv_isIso _ _ _ _ (rtExt P F A 1) (isIso_rtExt_one P F A)

/-- ★`𝒞` の側の代表元 `rtObj A 1 ⟶ rtObj A k`。 -/
noncomputable def kappaRep (A : C) (k : ℕ+) : rtObj P F A 1 ⟶ rtObj P F A k :=
  rtOneInv (P := P) (F := F) A ≫ rtExt P F A k

theorem kappaRep_degFr (A : C) (k : ℕ+) : P.degFr (kappaRep (P := P) (F := F) A k) = k := by
  have h1 : P.degFr (rtOneInv (P := P) (F := F) A) = 1 := degFr_of_isIso P _
  have h2 : P.degFr (rtExt P F A k) = k := rtExt_degFr P F A k
  show P.degFr (rtOneInv (P := P) (F := F) A ≫ rtExt P F A k) = k
  rw [P.degFr_comp, h1, h2, mul_one]

theorem kappaRep_frobType (A : C) (k : ℕ+) :
    IsFrobeniusType P (kappaRep (P := P) (F := F) A k) :=
  IsFrobeniusType.comp P F (isFrobeniusType_of_isIso P (rtOneInv (P := P) (F := F) A))
    (rtExt_frobType P F A k)

/-- ★★**標準の次数 `k` の Frobenius 型射 `⟨A,k⟩ ⟶ ⟨A,1⟩`**。 -/
noncomputable def pfKappa (A : C) (k : ℕ+) : (⟨A, k⟩ : PfRootObj P F) ⟶ ⟨A, 1⟩ :=
  toHomPf (F := F) (kappaRep (P := P) (F := F) A k)

theorem pfKappa_degFr (A : C) (k : ℕ+) :
    (pfRootPre P F).degFr (pfKappa (F := F) A k) = k := by
  show rootDeg (show HomRoot P F ⟨A, k⟩ ⟨A, 1⟩ from
    HomPf.mk (idxOne P F (rtObj P F A 1) (rtObj P F A k))
      (kappaRep (P := P) (F := F) A k)) = k
  rw [rootDeg_mk]
  exact kappaRep_degFr A k

theorem pfKappa_frobType (hfi : IsOfFrobeniusIsotropicType P) (A : C) (k : ℕ+) :
    IsFrobeniusType (pfRootPre P F) (pfKappa (F := F) A k) :=
  (isFrobeniusType_mk_iff (X := (⟨A, k⟩ : PfRootObj P F))
    (Y := (⟨A, 1⟩ : PfRootObj P F)) hfi
    (idxOne P F (rtObj P F A 1) (rtObj P F A k))
    (kappaRep (P := P) (F := F) A k)).mpr
    ⟨(kappaRep_frobType (P := P) (F := F) A k).1.2,
     (kappaRep_frobType (P := P) (F := F) A k).2⟩

/-- ★★**`κ` は mono** —— `𝒞^pf` の Frobenius 型射は mono(`pfRoot_frobTypeMono`)。 -/
theorem pfKappa_mono (hfi : IsOfFrobeniusIsotropicType P) (A : C) (k : ℕ+) :
    Mono (pfKappa (F := F) A k) :=
  pfRoot_frobTypeMono hfi _ (pfKappa_frobType hfi A k)

/-! ## ★2. `𝒞 → 𝒞^pf` の像の辞書(根 1 の形で) -/

theorem toRootHom_degFr' {A B : C} (φ : A ⟶ B) :
    (pfRootPre P F).degFr (toRootHom (F := F) φ) = P.degFr φ := toPfRoot_degFr φ

theorem toRootHom_frobType (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (φ : A ⟶ B)
    (h : IsFrobeniusType P φ) :
    IsFrobeniusType (pfRootPre P F) (toRootHom (F := F) φ) :=
  toPfRoot_isFrobeniusType hfi φ h

/-! ## ★3. `⟨A,1⟩ ≅ ⟨A₁,k⟩` —— 方程式で特徴づける -/

/-- ★★★★**持ち上げの存在** —— `e ≫ κ = [α]`。

★中身は `pfRoot_frob_div`(共通の終域への Frobenius 型射による割り算、在庫)。 -/
theorem exists_rootLift (hfi : IsOfFrobeniusIsotropicType P) {A A₁ : C} (α : A ⟶ A₁)
    (k : ℕ+) (hk : P.degFr α = k) :
    ∃ e : (⟨A, 1⟩ : PfRootObj P F) ⟶ ⟨A₁, k⟩,
      e ≫ pfKappa (F := F) A₁ k = toRootHom (F := F) α := by
  refine pfRoot_frob_div (F := F) (n := k) (m := 1) hfi (toRootHom (F := F) α)
    (pfKappa (F := F) A₁ k) ?_ (pfKappa_frobType hfi A₁ k) (pfKappa_degFr A₁ k)
  rw [toRootHom_degFr', hk, mul_one]

/-- ★★★★**持ち上げは同型** —— `κ` と `[α]` がどちらも mono だから。 -/
theorem isIso_rootLift (hfi : IsOfFrobeniusIsotropicType P) {A A₁ : C} (α : A ⟶ A₁)
    (hα : IsFrobeniusType P α) (k : ℕ+) (hk : P.degFr α = k)
    (e : (⟨A, 1⟩ : PfRootObj P F) ⟶ ⟨A₁, k⟩)
    (he : e ≫ pfKappa (F := F) A₁ k = toRootHom (F := F) α) : IsIso e := by
  obtain ⟨e', he'⟩ := pfRoot_frob_div (F := F) (n := k) (m := 1) hfi
    (pfKappa (F := F) A₁ k) (toRootHom (F := F) α)
    (by rw [pfKappa_degFr, mul_one])
    (toRootHom_frobType hfi α hα) (by rw [toRootHom_degFr', hk])
  haveI hm1 : Mono (toRootHom (F := F) α) :=
    pfRoot_frobTypeMono hfi _ (toRootHom_frobType hfi α hα)
  haveI hm2 : Mono (pfKappa (F := F) A₁ k) := pfKappa_mono hfi A₁ k
  refine ⟨e', ?_, ?_⟩
  · refine (cancel_mono (toRootHom (F := F) α)).mp ?_
    rw [Category.assoc, he', he, Category.id_comp]
  · refine (cancel_mono (pfKappa (F := F) A₁ k)).mp ?_
    rw [Category.assoc, he, he', Category.id_comp]

/-- ★★★**`⟨A,1⟩ ⟶ ⟨A₁,k⟩` の標準の射**。 -/
noncomputable def rootLift (hfi : IsOfFrobeniusIsotropicType P) {A A₁ : C} (α : A ⟶ A₁)
    (k : ℕ+) (hk : P.degFr α = k) : (⟨A, 1⟩ : PfRootObj P F) ⟶ ⟨A₁, k⟩ :=
  (exists_rootLift hfi α k hk).choose

@[simp] theorem rootLift_spec (hfi : IsOfFrobeniusIsotropicType P) {A A₁ : C} (α : A ⟶ A₁)
    (k : ℕ+) (hk : P.degFr α = k) :
    rootLift (F := F) hfi α k hk ≫ pfKappa (F := F) A₁ k = toRootHom (F := F) α :=
  (exists_rootLift hfi α k hk).choose_spec

theorem rootLift_isIso (hfi : IsOfFrobeniusIsotropicType P) {A A₁ : C} (α : A ⟶ A₁)
    (hα : IsFrobeniusType P α) (k : ℕ+) (hk : P.degFr α = k) :
    IsIso (rootLift (F := F) hfi α k hk) :=
  isIso_rootLift hfi α hα k hk _ (rootLift_spec hfi α k hk)

/-! ## ★4. 根を上げる標準の射 -/

/-- ★★★**根を上げる標準の射の存在** —— `E ≫ κ′ = κ ≫ [a]`。 -/
theorem exists_rootStep (hfi : IsOfFrobeniusIsotropicType P) {A₁ A₁' : C} (a : A₁ ⟶ A₁')
    (ha : IsFrobeniusType P a) (k d k' : ℕ+) (hd : P.degFr a = d) (hk' : k' = k * d) :
    ∃ E : (⟨A₁, k⟩ : PfRootObj P F) ⟶ ⟨A₁', k'⟩,
      E ≫ pfKappa (F := F) A₁' k' = pfKappa (F := F) A₁ k ≫ toRootHom (F := F) a := by
  refine pfRoot_frob_div (F := F) (n := k') (m := 1) hfi
    (pfKappa (F := F) A₁ k ≫ toRootHom (F := F) a) (pfKappa (F := F) A₁' k') ?_
    (pfKappa_frobType hfi A₁' k') (pfKappa_degFr A₁' k')
  rw [(pfRootPre P F).degFr_comp, toRootHom_degFr', pfKappa_degFr, hd, hk', mul_one, mul_comm]

/-- ★★★**根を上げる標準の射は同型**。 -/
theorem isIso_rootStep (hfi : IsOfFrobeniusIsotropicType P) {A₁ A₁' : C} (a : A₁ ⟶ A₁')
    (ha : IsFrobeniusType P a) (k d k' : ℕ+) (hd : P.degFr a = d) (hk' : k' = k * d)
    (E : (⟨A₁, k⟩ : PfRootObj P F) ⟶ ⟨A₁', k'⟩)
    (hE : E ≫ pfKappa (F := F) A₁' k' = pfKappa (F := F) A₁ k ≫ toRootHom (F := F) a) :
    IsIso E := by
  have hcomp : IsFrobeniusType (pfRootPre P F)
      (pfKappa (F := F) A₁ k ≫ toRootHom (F := F) a) :=
    IsFrobeniusType.comp (pfRootPre P F) (pfRootCore hfi) (pfKappa_frobType hfi A₁ k)
      (toRootHom_frobType hfi a ha)
  have hdeg : (pfRootPre P F).degFr (pfKappa (F := F) A₁ k ≫ toRootHom (F := F) a) = k' := by
    rw [(pfRootPre P F).degFr_comp, toRootHom_degFr', pfKappa_degFr, hd, hk', mul_comm]
  obtain ⟨E', hE'⟩ := pfRoot_frob_div (F := F) (n := k') (m := 1) hfi
    (pfKappa (F := F) A₁' k') (pfKappa (F := F) A₁ k ≫ toRootHom (F := F) a)
    (by rw [pfKappa_degFr, mul_one]) hcomp hdeg
  haveI hm1 : Mono (pfKappa (F := F) A₁ k ≫ toRootHom (F := F) a) :=
    pfRoot_frobTypeMono hfi _ hcomp
  haveI hm2 : Mono (pfKappa (F := F) A₁' k') := pfKappa_mono hfi A₁' k'
  refine ⟨E', ?_, ?_⟩
  · refine (cancel_mono (pfKappa (F := F) A₁ k ≫ toRootHom (F := F) a)).mp ?_
    rw [Category.assoc, hE', hE, Category.id_comp]
  · refine (cancel_mono (pfKappa (F := F) A₁' k')).mp ?_
    rw [Category.assoc, hE, hE', Category.id_comp]

/-- ★★★**根を上げる標準の射**。 -/
noncomputable def rootStep (hfi : IsOfFrobeniusIsotropicType P) {A₁ A₁' : C} (a : A₁ ⟶ A₁')
    (ha : IsFrobeniusType P a) (k d k' : ℕ+) (hd : P.degFr a = d) (hk' : k' = k * d) :
    (⟨A₁, k⟩ : PfRootObj P F) ⟶ ⟨A₁', k'⟩ :=
  (exists_rootStep (F := F) hfi a ha k d k' hd hk').choose

@[simp] theorem rootStep_spec (hfi : IsOfFrobeniusIsotropicType P) {A₁ A₁' : C} (a : A₁ ⟶ A₁')
    (ha : IsFrobeniusType P a) (k d k' : ℕ+) (hd : P.degFr a = d) (hk' : k' = k * d) :
    rootStep (F := F) hfi a ha k d k' hd hk' ≫ pfKappa (F := F) A₁' k'
      = pfKappa (F := F) A₁ k ≫ toRootHom (F := F) a :=
  (exists_rootStep (F := F) hfi a ha k d k' hd hk').choose_spec

theorem rootStep_isIso (hfi : IsOfFrobeniusIsotropicType P) {A₁ A₁' : C} (a : A₁ ⟶ A₁')
    (ha : IsFrobeniusType P a) (k d k' : ℕ+) (hd : P.degFr a = d) (hk' : k' = k * d) :
    IsIso (rootStep (F := F) hfi a ha k d k' hd hk') :=
  isIso_rootStep hfi a ha k d k' hd hk' _ (rootStep_spec hfi a ha k d k' hd hk')

/-! ## ★5. 整合性 -/

/-- ★★★★★★**整合性** —— `rootLift (α ≫ a) = rootLift α ≫ rootStep a`。

★★これが「選択された同型」の問題を消す 1 本である ——
`κ` が mono なので、`κ` と合成した式が一致すれば射そのものが一致する。 -/
theorem rootLift_comp (hfi : IsOfFrobeniusIsotropicType P) {A A₁ A₁' : C}
    (α : A ⟶ A₁) (a : A₁ ⟶ A₁') (ha : IsFrobeniusType P a)
    (k d k' : ℕ+) (hk : P.degFr α = k) (hd : P.degFr a = d) (hk' : k' = k * d)
    (hka : P.degFr (α ≫ a) = k') :
    rootLift (F := F) hfi (α ≫ a) k' hka
      = rootLift (F := F) hfi α k hk ≫ rootStep (F := F) hfi a ha k d k' hd hk' := by
  haveI hm : Mono (pfKappa (F := F) A₁' k') := pfKappa_mono hfi A₁' k'
  refine (cancel_mono (pfKappa (F := F) A₁' k')).mp ?_
  rw [rootLift_spec, Category.assoc, rootStep_spec, ← Category.assoc, rootLift_spec]
  exact toRootHom_comp α a

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Proposition 5.5, (ii)` で選択された同型の整合性を消す 1 本。 -/
def rootLift_comp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 𝒞^pf の根の標準射と整合性",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
