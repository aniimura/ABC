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

/-! ## ★6. 根を上げた射 —— こちらも方程式で特徴づける

★★`lamHom`(在庫)は `⟨A,k⟩ ⟶ ⟨B,k⟩` を与えるが、
**`κ` との自然性を直接計算しようとすると `compRoot` の展開が要る**。
★そこで `rootMap` を**方程式 `g ≫ κ_B = κ_A ≫ [f]` で定義する** ——
存在は `pfRoot_frob_div`、一意性は `κ_B` が mono であることから出る。
★★これで自然性も関手性も**両辺に `κ` を合成して比べるだけ**になる。 -/

theorem toRootHom_baseIso {A B : C} (f : A ⟶ B) (hf : IsBaseIsomorphism P f) :
    IsBaseIsomorphism (pfRootPre P F) (toRootHom (F := F) f) := by
  have hb1 : IsBaseIsomorphism P (rtOneInv (P := P) (F := F) A) :=
    isBaseIsomorphism_of_isIso P _
  have hb2 : IsBaseIsomorphism P (f ≫ rtExt P F B 1) :=
    isBaseIsomorphism_comp P hf (rtExt_frobType P F B 1).2
  exact (isBaseIsomorphism_mk_iff (X := (⟨A, 1⟩ : PfRootObj P F))
    (Y := (⟨B, 1⟩ : PfRootObj P F))
    (idxOne P F (rtObj P F A 1) (rtObj P F B 1))
    (rtOneInv (P := P) (F := F) A ≫ f ≫ rtExt P F B 1)).mpr
    (isBaseIsomorphism_comp P hb1 hb2)

/-- ★★★**根を上げた射の存在** —— `g ≫ κ_B = κ_A ≫ [f]`。 -/
theorem exists_rootMap (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (f : A ⟶ B) (k : ℕ+) :
    ∃ g : (⟨A, k⟩ : PfRootObj P F) ⟶ ⟨B, k⟩,
      g ≫ pfKappa (F := F) B k = pfKappa (F := F) A k ≫ toRootHom (F := F) f := by
  refine pfRoot_frob_div (F := F) (n := k) (m := P.degFr f) hfi
    (pfKappa (F := F) A k ≫ toRootHom (F := F) f) (pfKappa (F := F) B k) ?_
    (pfKappa_frobType hfi B k) (pfKappa_degFr B k)
  rw [(pfRootPre P F).degFr_comp, toRootHom_degFr', pfKappa_degFr, mul_comm]

/-- ★★**根を上げた射**(方程式つき)。 -/
noncomputable def rootMap (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (f : A ⟶ B) (k : ℕ+) :
    (⟨A, k⟩ : PfRootObj P F) ⟶ ⟨B, k⟩ :=
  (exists_rootMap (F := F) hfi f k).choose

@[simp] theorem rootMap_spec (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (f : A ⟶ B)
    (k : ℕ+) :
    rootMap (F := F) hfi f k ≫ pfKappa (F := F) B k
      = pfKappa (F := F) A k ≫ toRootHom (F := F) f :=
  (exists_rootMap (F := F) hfi f k).choose_spec

/-- ★★**`⟨B,k⟩` への射は `κ_B` と合成すれば決まる**(`κ` が mono)。 -/
theorem rootMap_ext (hfi : IsOfFrobeniusIsotropicType P) {X : PfRootObj P F} {B : C} {k : ℕ+}
    (g g' : X ⟶ (⟨B, k⟩ : PfRootObj P F))
    (h : g ≫ pfKappa (F := F) B k = g' ≫ pfKappa (F := F) B k) : g = g' := by
  haveI := pfKappa_mono (F := F) hfi B k
  exact (cancel_mono (pfKappa (F := F) B k)).mp h

theorem rootMap_degFr (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (f : A ⟶ B) (k : ℕ+) :
    (pfRootPre P F).degFr (rootMap (F := F) hfi f k) = P.degFr f := by
  have h := congrArg (pfRootPre P F).degFr (rootMap_spec (F := F) hfi f k)
  rw [(pfRootPre P F).degFr_comp, (pfRootPre P F).degFr_comp, pfKappa_degFr, pfKappa_degFr,
    toRootHom_degFr'] at h
  exact mul_left_cancel (a := k) (h.trans (mul_comm _ _))

/-- ★★**関手性**。 -/
theorem rootMap_comp (hfi : IsOfFrobeniusIsotropicType P) {A B E : C} (f : A ⟶ B) (g : B ⟶ E)
    (k : ℕ+) :
    rootMap (F := F) hfi (f ≫ g) k
      = rootMap (F := F) hfi f k ≫ rootMap (F := F) hfi g k := by
  refine rootMap_ext hfi _ _ ?_
  rw [rootMap_spec, Category.assoc, rootMap_spec, ← Category.assoc, rootMap_spec,
    Category.assoc]
  exact congrArg (fun t => pfKappa (F := F) A k ≫ t) (toRootHom_comp f g)

theorem rootMap_baseIso (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (f : A ⟶ B) (k : ℕ+)
    (hf : IsBaseIsomorphism P f) :
    IsBaseIsomorphism (pfRootPre P F) (rootMap (F := F) hfi f k) := by
  have h := congrArg (pfRootPre P F).Base (rootMap_spec (F := F) hfi f k)
  rw [(pfRootPre P F).Base_comp, (pfRootPre P F).Base_comp] at h
  haveI h1 : IsIso ((pfRootPre P F).Base (pfKappa (F := F) B k)) :=
    (pfKappa_frobType hfi B k).2
  haveI h2 : IsIso ((pfRootPre P F).Base (pfKappa (F := F) A k)) :=
    (pfKappa_frobType hfi A k).2
  haveI h3 : IsIso ((pfRootPre P F).Base (toRootHom (F := F) f)) := toRootHom_baseIso f hf
  haveI h4 : IsIso ((pfRootPre P F).Base (rootMap (F := F) hfi f k)
      ≫ (pfRootPre P F).Base (pfKappa (F := F) B k)) := h ▸ inferInstance
  exact IsIso.of_isIso_comp_right _ ((pfRootPre P F).Base (pfKappa (F := F) B k))

/-- ★★**pre-step は保たれる**。 -/
theorem rootMap_preStep (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (f : A ⟶ B) (k : ℕ+)
    (hf : IsPreStep P f) :
    IsPreStep (pfRootPre P F) (rootMap (F := F) hfi f k) :=
  ⟨show (pfRootPre P F).degFr (rootMap (F := F) hfi f k) = 1 from
    (rootMap_degFr hfi f k).trans hf.1,
   rootMap_baseIso hfi f k hf.2⟩

/-- ★★★★★**`rootMap` と `rootStep` の自然性の四角形**。

★`𝒞` の可換な四角形 `z ≫ a = β ≫ α`(`a`, `β` は同次数の Frobenius 型)を
`𝒞^pf` の根つきの四角形へ移す。★★中身は**両辺に `κ` を合成して比べるだけ**。 -/
theorem rootMap_rootStep_sq (hfi : IsOfFrobeniusIsotropicType P) {X W₁ Y W₁' : C}
    (z : X ⟶ W₁) (a : W₁ ⟶ W₁') (ha : IsFrobeniusType P a)
    (β : X ⟶ Y) (hβ : IsFrobeniusType P β) (α : Y ⟶ W₁')
    (k d k' : ℕ+) (hda : P.degFr a = d) (hdβ : P.degFr β = d) (hk' : k' = k * d)
    (hsq : z ≫ a = β ≫ α) :
    rootMap (F := F) hfi z k ≫ rootStep (F := F) hfi a ha k d k' hda hk'
      = rootStep (F := F) hfi β hβ k d k' hdβ hk' ≫ rootMap (F := F) hfi α k' := by
  refine rootMap_ext hfi _ _ ?_
  rw [Category.assoc, rootStep_spec, ← Category.assoc, rootMap_spec, Category.assoc,
    Category.assoc, rootMap_spec]
  have hR : rootStep (F := F) hfi β hβ k d k' hdβ hk'
        ≫ pfKappa (F := F) Y k' ≫ toRootHom (F := F) α
      = (pfKappa (F := F) X k ≫ toRootHom (F := F) β) ≫ toRootHom (F := F) α := by
    rw [← Category.assoc, rootStep_spec]
  rw [hR, Category.assoc]
  exact congrArg (fun t => pfKappa (F := F) X k ≫ t)
    ((toRootHom_comp z a).symm.trans ((congrArg toRootHom hsq).trans (toRootHom_comp β α)))

/-- ★★★★locator —— `Proposition 5.5, (ii)` の根つきの自然性。 -/
def rootMap_rootStep_sq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 根を上げた射の自然性",
    sectionId := "frdi-prop-5-5" }

/-! ## ★7. `Hom^pf` の代表元を「標準の `k` 乗根の添字」に揃える

★★`Proposition 5.5, (ii)` の**逆写像**を作るときに要る ——
左辺の元は「`𝒞^pf` の co-angular pre-step `ε`」と「`𝒞^pf` の射 `φ`」の組であり、
どちらも `Hom^pf` の元なので、**共通の根に揃える**と
`𝒞` のスパン `A^{(k)} ← X → B^{(k)}` になる。

★添字 `Z` の次数を `k` とすると `Definition 1.3, (ii)` の本質的一意性
(`frobDegUniq`)が `Z ≅ idxPow U V k` を与えるので、
**どの元も `idxPow` の形の添字で代表できる**。 -/

/-- ★★★**`Hom^pf` の元は「標準の `k` 乗根の添字」で代表できる**。 -/
theorem HomPf.exists_rep_pow {U V : C} (z : HomPf P F U V) :
    ∃ (k : ℕ+) (ψ : rtObj P F U k ⟶ rtObj P F V k),
      HomPf.mk (idxPow (F := F) U V k) ψ = z := by
  obtain ⟨Z, φ, hφ⟩ := HomPf.exists_rep z
  obtain ⟨θ₁, hθ₁, hθ₁e⟩ := F.frobDegUniq U Z.right.obj.1
    (rtObj P F U (P.degFr Z.hom.hom.1))
    Z.hom.hom.1 (rtExt P F U (P.degFr Z.hom.hom.1)) Z.hom.property.1
    (rtExt_frobType P F U _) (rtExt_degFr P F U _).symm
  obtain ⟨θ₂, hθ₂, hθ₂e⟩ := F.frobDegUniq V Z.right.obj.2
    (rtObj P F V (P.degFr Z.hom.hom.1))
    Z.hom.hom.2 (rtExt P F V (P.degFr Z.hom.hom.1)) Z.hom.property.2.1
    (rtExt_frobType P F V _)
    (Z.hom.property.2.2.symm.trans (rtExt_degFr P F V _).symm)
  haveI := hθ₁
  haveI := hθ₂
  have hd1 : P.degFr θ₁ = 1 := degFr_of_isIso P θ₁
  have hd2 : P.degFr θ₂ = 1 := degFr_of_isIso P θ₂
  refine ⟨P.degFr Z.hom.hom.1, idxTransport P F (Under.homMk
    (show Z.right ⟶ (idxPow (F := F) U V (P.degFr Z.hom.hom.1)).right from
      ⟨(θ₁, θ₂), isFrobeniusType_of_isIso P θ₁, isFrobeniusType_of_isIso P θ₂,
        hd1.trans hd2.symm⟩)
    (WideSubcategory.hom_ext _ (Prod.ext hθ₁e hθ₂e))) φ, ?_⟩
  rw [← hφ]
  exact HomPf.mk_map _ φ

/-- ★添字 `idxPow U V k ⟶ idxPow U V t`(`t = e * k`)。 -/
noncomputable def idxPowLift (U V : C) {k e t : ℕ+} (ht : t = e * k) :
    idxPow (F := F) U V k ⟶ idxPow (F := F) U V t :=
  Under.homMk
    (show (idxPow (F := F) U V k).right ⟶ (idxPow (F := F) U V t).right from
      ⟨(rtLift P F U ht, rtLift P F V ht), rtLift_frobType P F U ht,
        rtLift_frobType P F V ht,
        (rtLift_degFr P F U ht).trans (rtLift_degFr P F V ht).symm⟩)
    (WideSubcategory.hom_ext _ (Prod.ext (rtLift_ext P F U ht) (rtLift_ext P F V ht)))

/-- ★★**代表元は根を上げても同じ元を与える**。 -/
theorem HomPf.mk_pow_lift {U V : C} {k e t : ℕ+} (ht : t = e * k)
    (ψ : rtObj P F U k ⟶ rtObj P F V k) :
    HomPf.mk (idxPow (F := F) U V t) (idxTransport P F (idxPowLift (F := F) U V ht) ψ)
      = HomPf.mk (idxPow (F := F) U V k) ψ :=
  HomPf.mk_map (idxPowLift (F := F) U V ht) ψ

/-- ★★★**同じ始域から出る 2 つの `Hom^pf` の元は共通の根で代表できる**。

★★これが逆写像の第 1 歩である。 -/
theorem HomPf.exists_rep_pow_pair {U X Y : C} (z₁ : HomPf P F U X) (z₂ : HomPf P F U Y) :
    ∃ (k : ℕ+) (ψ₁ : rtObj P F U k ⟶ rtObj P F X k) (ψ₂ : rtObj P F U k ⟶ rtObj P F Y k),
      HomPf.mk (idxPow (F := F) U X k) ψ₁ = z₁ ∧
      HomPf.mk (idxPow (F := F) U Y k) ψ₂ = z₂ := by
  obtain ⟨k₁, ψ₁, h₁⟩ := HomPf.exists_rep_pow z₁
  obtain ⟨k₂, ψ₂, h₂⟩ := HomPf.exists_rep_pow z₂
  refine ⟨k₂ * k₁,
    idxTransport P F (idxPowLift (F := F) U X (e := k₂) (t := k₂ * k₁) rfl) ψ₁,
    idxTransport P F (idxPowLift (F := F) U Y (e := k₁) (t := k₂ * k₁)
      (mul_comm k₂ k₁)) ψ₂, ?_, ?_⟩
  · exact (HomPf.mk_pow_lift (F := F) rfl ψ₁).trans h₁
  · exact (HomPf.mk_pow_lift (F := F) (mul_comm k₂ k₁) ψ₂).trans h₂

/-- ★★★locator —— `Proposition 5.5, (ii)` の逆写像の第 1 歩(共通の根)。 -/
def HomPf.exists_rep_pow_pair.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Hom^pf の 2 元を共通の根で代表する",
    sectionId := "frdi-prop-5-5" }

/-! ## ★8. `κ` の一般形 —— 始域も終域も根つきの場合

★★`Proposition 5.5, (ii)` の**逆向き**では、始域が `⟨A″, n⟩` の形になる。
★`A″` の `n` 乗根は `𝒞` に無いので、根 1 を終域に持つ `pfKappa` だけでは足りない。
★★そこで **`κ` を `⟨A, k·r⟩ ⟶ ⟨A, r⟩` へ一般化する** ——
代表元は `rtLift`(在庫、次数 `k` の Frobenius 型)そのものである。
★これで `rootLift` / `rootStep` / `rootMap` の 3 本がそのまま一般化する。 -/

/-- ★sharp な単系では `k • x = 0` から `x = 0`。 -/
theorem eq_zero_of_nsmul_eq_zero_sharp {M : Type*} [AddCommMonoid M] (hs : IsSharp M)
    {x : M} {k : ℕ} (hk : 0 < k) (h : k • x = 0) : x = 0 := by
  refine hs x ⟨⟨x, (k - 1) • x, ?_, ?_⟩, rfl⟩
  · have h1 : x + (k - 1) • x = (1 + (k - 1)) • x := by rw [add_nsmul, one_nsmul]
    rw [h1, show 1 + (k - 1) = k from by omega]
    exact h
  · have h1 : (k - 1) • x + x = ((k - 1) + 1) • x := by rw [add_nsmul, one_nsmul]
    rw [h1, show (k - 1) + 1 = k from by omega]
    exact h

/-- ★★**`κ` の一般形** `⟨A, t⟩ ⟶ ⟨A, r⟩`(`t = k * r`)。 -/
noncomputable def pfKappaGen (A : C) {k r t : ℕ+} (ht : t = k * r) :
    (⟨A, t⟩ : PfRootObj P F) ⟶ ⟨A, r⟩ :=
  toHomPf (F := F) (rtLift P F A ht)

theorem pfKappaGen_degFr (A : C) {k r t : ℕ+} (ht : t = k * r) :
    (pfRootPre P F).degFr (pfKappaGen (F := F) A ht) = k := by
  show rootDeg (show HomRoot P F ⟨A, t⟩ ⟨A, r⟩ from
    HomPf.mk (idxOne P F (rtObj P F A r) (rtObj P F A t)) (rtLift P F A ht)) = k
  rw [rootDeg_mk]
  exact rtLift_degFr P F A ht

theorem pfKappaGen_frobType (hfi : IsOfFrobeniusIsotropicType P) (A : C) {k r t : ℕ+}
    (ht : t = k * r) : IsFrobeniusType (pfRootPre P F) (pfKappaGen (F := F) A ht) :=
  (isFrobeniusType_mk_iff (X := (⟨A, t⟩ : PfRootObj P F))
    (Y := (⟨A, r⟩ : PfRootObj P F)) hfi
    (idxOne P F (rtObj P F A r) (rtObj P F A t)) (rtLift P F A ht)).mpr
    ⟨(rtLift_frobType P F A ht).1.2, (rtLift_frobType P F A ht).2⟩

theorem pfKappaGen_mono (hfi : IsOfFrobeniusIsotropicType P) (A : C) {k r t : ℕ+}
    (ht : t = k * r) : Mono (pfKappaGen (F := F) A ht) :=
  pfRoot_frobTypeMono hfi _ (pfKappaGen_frobType hfi A ht)

/-- ★★**`rootMap` は Frobenius 型射に対して mono**。

★`rootMap α r ≫ κ = κ ≫ [α]` の右辺が Frobenius 型(ゆえに mono)だから。 -/
theorem rootMap_mono (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (α : A ⟶ B)
    (hα : IsFrobeniusType P α) (r : ℕ+) : Mono (rootMap (F := F) hfi α r) := by
  have hcomp : IsFrobeniusType (pfRootPre P F)
      (pfKappa (F := F) A r ≫ toRootHom (F := F) α) :=
    IsFrobeniusType.comp (pfRootPre P F) (pfRootCore hfi) (pfKappa_frobType hfi A r)
      (toRootHom_frobType hfi α hα)
  haveI hm : Mono (pfKappa (F := F) A r ≫ toRootHom (F := F) α) :=
    pfRoot_frobTypeMono hfi _ hcomp
  haveI hm2 : Mono (rootMap (F := F) hfi α r ≫ pfKappa (F := F) B r) := by
    rw [rootMap_spec]; exact hm
  exact mono_of_mono (rootMap (F := F) hfi α r) (pfKappa (F := F) B r)

/-- ★★**`rootMap` は等長性を保つ** —— `Div_comp` と sharp から。 -/
theorem rootMap_isometric (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (α : A ⟶ B)
    (hα : IsFrobeniusType P α) (r : ℕ+) :
    IsIsometric (pfRootPre P F) (rootMap (F := F) hfi α r) := by
  have h := congrArg (pfRootPre P F).Div (rootMap_spec (F := F) hfi α r)
  rw [(pfRootPre P F).Div_comp, (pfRootPre P F).Div_comp] at h
  have h1 : (pfRootPre P F).Div (pfKappa (F := F) B r) = 0 :=
    (pfKappa_frobType hfi B r).1.2
  have h2 : (pfRootPre P F).Div (pfKappa (F := F) A r) = 0 :=
    (pfKappa_frobType hfi A r).1.2
  have h3 : (pfRootPre P F).Div (toRootHom (F := F) α) = 0 :=
    (toRootHom_frobType hfi α hα).1.2
  rw [h1, h2, h3, map_zero, map_zero, smul_zero, zero_add, add_zero] at h
  rw [pfKappa_degFr] at h
  exact eq_zero_of_nsmul_eq_zero_sharp ((pfRootPre P F).divisorial _).2 (k := (r : ℕ)) r.2 h

/-- ★★**`rootMap` は Frobenius 型を保つ**。 -/
theorem rootMap_frobType (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (α : A ⟶ B)
    (hα : IsFrobeniusType P α) (r : ℕ+) :
    IsFrobeniusType (pfRootPre P F) (rootMap (F := F) hfi α r) :=
  ⟨⟨pfRoot_isCoAngular hfi _, rootMap_isometric hfi α hα r⟩, rootMap_baseIso hfi α r hα.2⟩

/-- ★★★★**根つきの持ち上げの存在** —— `e ≫ κ_gen = rootMap α r`。 -/
theorem exists_rootLiftGen (hfi : IsOfFrobeniusIsotropicType P) {A A₁ : C} (α : A ⟶ A₁)
    {k r t : ℕ+} (hk : P.degFr α = k) (ht : t = k * r) :
    ∃ e : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A₁, t⟩,
      e ≫ pfKappaGen (F := F) A₁ ht = rootMap (F := F) hfi α r := by
  refine pfRoot_frob_div (F := F) (n := k) (m := 1) hfi (rootMap (F := F) hfi α r)
    (pfKappaGen (F := F) A₁ ht) ?_ (pfKappaGen_frobType hfi A₁ ht) (pfKappaGen_degFr A₁ ht)
  rw [rootMap_degFr, hk, mul_one]

/-- ★★★★**根つきの持ち上げは同型**。 -/
theorem isIso_rootLiftGen (hfi : IsOfFrobeniusIsotropicType P) {A A₁ : C} (α : A ⟶ A₁)
    (hα : IsFrobeniusType P α) {k r t : ℕ+} (hk : P.degFr α = k) (ht : t = k * r)
    (e : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A₁, t⟩)
    (he : e ≫ pfKappaGen (F := F) A₁ ht = rootMap (F := F) hfi α r) : IsIso e := by
  obtain ⟨e', he'⟩ := pfRoot_frob_div (F := F) (n := k) (m := 1) hfi
    (pfKappaGen (F := F) A₁ ht) (rootMap (F := F) hfi α r)
    (by rw [pfKappaGen_degFr, mul_one])
    (rootMap_frobType hfi α hα r) (by rw [rootMap_degFr, hk])
  haveI hm1 : Mono (rootMap (F := F) hfi α r) := rootMap_mono hfi α hα r
  haveI hm2 : Mono (pfKappaGen (F := F) A₁ ht) := pfKappaGen_mono hfi A₁ ht
  refine ⟨e', ?_, ?_⟩
  · refine (cancel_mono (rootMap (F := F) hfi α r)).mp ?_
    rw [Category.assoc, he', he, Category.id_comp]
  · refine (cancel_mono (pfKappaGen (F := F) A₁ ht)).mp ?_
    rw [Category.assoc, he, he', Category.id_comp]

/-- ★★★**根つきの標準の持ち上げ** `⟨A,r⟩ ≅ ⟨A₁, k·r⟩`。 -/
noncomputable def rootLiftGen (hfi : IsOfFrobeniusIsotropicType P) {A A₁ : C} (α : A ⟶ A₁)
    {k r t : ℕ+} (hk : P.degFr α = k) (ht : t = k * r) :
    (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A₁, t⟩ :=
  (exists_rootLiftGen hfi α hk ht).choose

@[simp] theorem rootLiftGen_spec (hfi : IsOfFrobeniusIsotropicType P) {A A₁ : C} (α : A ⟶ A₁)
    {k r t : ℕ+} (hk : P.degFr α = k) (ht : t = k * r) :
    rootLiftGen (F := F) hfi α hk ht ≫ pfKappaGen (F := F) A₁ ht
      = rootMap (F := F) hfi α r :=
  (exists_rootLiftGen hfi α hk ht).choose_spec

theorem rootLiftGen_isIso (hfi : IsOfFrobeniusIsotropicType P) {A A₁ : C} (α : A ⟶ A₁)
    (hα : IsFrobeniusType P α) {k r t : ℕ+} (hk : P.degFr α = k) (ht : t = k * r) :
    IsIso (rootLiftGen (F := F) hfi α hk ht) :=
  isIso_rootLiftGen hfi α hα hk ht _ (rootLiftGen_spec hfi α hk ht)

/-- ★★★★locator —— `Proposition 5.5, (ii)` の逆向きに要る `κ` の一般形。 -/
def rootLiftGen.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 根つきの κ と持ち上げ",
    sectionId := "frdi-prop-5-5" }

/-! ## ★9. 押し出しの添字 —— `rootIso` の計算則に直結する形

★★★**実務上の要点**: `idxPow U V k`(`idxMk (rtExt U k) (rtExt V k)`)と
「押し出し」`(pushIdx (rtExt U k) … (rtExt V k) …).obj (idxOne …)` は
**構造射が `≫ 𝟙` の分だけ違う**。★`≫ 𝟙` は一般の圏では `rfl` ではないので、
2 つの添字対象は**定義的には等しくない**。

★★`WideSubcategory` の恒等射は `.hom` 越しに簡約されないため、
「両者を結ぶ自明な添字の射」を作る道は詰まる(実測)。
★★★**押し出しの形をそのまま使うのが正解**である ——
`rootIso_hom_mk`(在庫)が**そのまま**当たる。 -/

/-- ★`idxPow` を「押し出し」の形で書いたもの(`rootIso` の計算則に直結する)。 -/
noncomputable def idxPow' (U V : C) (k : ℕ+) : IdxPf P F U V :=
  (pushIdx (F := F) (rtExt P F U k) (rtExt_frobType P F U k)
    (rtExt P F V k) (rtExt_frobType P F V k)
    (by rw [rtExt_degFr, rtExt_degFr])).obj (idxOne P F (rtObj P F U k) (rtObj P F V k))

/-- ★★★**辞書** —— 押し出しの添字での代表元は `rootIso` で `toHomPf` に戻る。

★★これが `Proposition 5.5, (ii)` の**逆向き**に要る 1 本である。 -/
theorem HomPf.mk_idxPow' (U V : C) (k : ℕ+) (ψ : rtObj P F U k ⟶ rtObj P F V k) :
    HomPf.mk (idxPow' (F := F) U V k) ψ
      = (rootIso (F := F) (rtExt P F U k) (rtExt_frobType P F U k)
          (rtExt P F V k) (rtExt_frobType P F V k)
          (by rw [rtExt_degFr, rtExt_degFr])).hom (toHomPf (F := F) ψ) :=
  (rootIso_hom_mk (F := F) (rtExt P F U k) (rtExt_frobType P F U k)
    (rtExt P F V k) (rtExt_frobType P F V k)
    (by rw [rtExt_degFr, rtExt_degFr])
    (idxOne P F (rtObj P F U k) (rtObj P F V k)) ψ).symm

/-- ★★★**`Hom^pf` の元は「押し出しの添字」で代表できる**。 -/
theorem HomPf.exists_rep_pow' {U V : C} (z : HomPf P F U V) :
    ∃ (k : ℕ+) (ψ : rtObj P F U k ⟶ rtObj P F V k),
      HomPf.mk (idxPow' (F := F) U V k) ψ = z := by
  obtain ⟨Z, φ, hφ⟩ := HomPf.exists_rep z
  obtain ⟨θ₁, hθ₁, hθ₁e⟩ := F.frobDegUniq U Z.right.obj.1
    (rtObj P F U (P.degFr Z.hom.hom.1))
    Z.hom.hom.1 (rtExt P F U (P.degFr Z.hom.hom.1)) Z.hom.property.1
    (rtExt_frobType P F U _) (rtExt_degFr P F U _).symm
  obtain ⟨θ₂, hθ₂, hθ₂e⟩ := F.frobDegUniq V Z.right.obj.2
    (rtObj P F V (P.degFr Z.hom.hom.1))
    Z.hom.hom.2 (rtExt P F V (P.degFr Z.hom.hom.1)) Z.hom.property.2.1
    (rtExt_frobType P F V _)
    (Z.hom.property.2.2.symm.trans (rtExt_degFr P F V _).symm)
  haveI := hθ₁
  haveI := hθ₂
  have hd1 : P.degFr θ₁ = 1 := degFr_of_isIso P θ₁
  have hd2 : P.degFr θ₂ = 1 := degFr_of_isIso P θ₂
  refine ⟨P.degFr Z.hom.hom.1, idxTransport P F (Under.homMk
    (show Z.right ⟶ (idxPow' (F := F) U V (P.degFr Z.hom.hom.1)).right from
      ⟨(θ₁, θ₂), isFrobeniusType_of_isIso P θ₁, isFrobeniusType_of_isIso P θ₂,
        hd1.trans hd2.symm⟩)
    (WideSubcategory.hom_ext _
      (Prod.ext (hθ₁e.trans (Category.comp_id _).symm)
        (hθ₂e.trans (Category.comp_id _).symm)))) φ, ?_⟩
  rw [← hφ]
  exact HomPf.mk_map _ φ

/-- ★添字 `idxPow' U V k ⟶ idxPow' U V t`(`t = e * k`)。 -/
noncomputable def idxPowLift' (U V : C) {k e t : ℕ+} (ht : t = e * k) :
    idxPow' (F := F) U V k ⟶ idxPow' (F := F) U V t :=
  Under.homMk
    (show (idxPow' (F := F) U V k).right ⟶ (idxPow' (F := F) U V t).right from
      ⟨(rtLift P F U ht, rtLift P F V ht), rtLift_frobType P F U ht,
        rtLift_frobType P F V ht,
        (rtLift_degFr P F U ht).trans (rtLift_degFr P F V ht).symm⟩)
    (WideSubcategory.hom_ext _
      (Prod.ext
        (show (rtExt P F U k ≫ 𝟙 (rtObj P F U k)) ≫ rtLift P F U ht
            = rtExt P F U t ≫ 𝟙 (rtObj P F U t) by
          rw [Category.comp_id, Category.comp_id]; exact rtLift_ext P F U ht)
        (show (rtExt P F V k ≫ 𝟙 (rtObj P F V k)) ≫ rtLift P F V ht
            = rtExt P F V t ≫ 𝟙 (rtObj P F V t) by
          rw [Category.comp_id, Category.comp_id]; exact rtLift_ext P F V ht)))

theorem HomPf.mk_pow_lift' {U V : C} {k e t : ℕ+} (ht : t = e * k)
    (ψ : rtObj P F U k ⟶ rtObj P F V k) :
    HomPf.mk (idxPow' (F := F) U V t) (idxTransport P F (idxPowLift' (F := F) U V ht) ψ)
      = HomPf.mk (idxPow' (F := F) U V k) ψ :=
  HomPf.mk_map (idxPowLift' (F := F) U V ht) ψ

/-- ★★★**同じ始域から出る 2 つの `Hom^pf` の元は共通の根で代表できる**(押し出し版)。 -/
theorem HomPf.exists_rep_pow_pair' {U X Y : C} (z₁ : HomPf P F U X) (z₂ : HomPf P F U Y) :
    ∃ (k : ℕ+) (ψ₁ : rtObj P F U k ⟶ rtObj P F X k) (ψ₂ : rtObj P F U k ⟶ rtObj P F Y k),
      HomPf.mk (idxPow' (F := F) U X k) ψ₁ = z₁ ∧
      HomPf.mk (idxPow' (F := F) U Y k) ψ₂ = z₂ := by
  obtain ⟨k₁, ψ₁, h₁⟩ := HomPf.exists_rep_pow' z₁
  obtain ⟨k₂, ψ₂, h₂⟩ := HomPf.exists_rep_pow' z₂
  refine ⟨k₂ * k₁,
    idxTransport P F (idxPowLift' (F := F) U X (e := k₂) (t := k₂ * k₁) rfl) ψ₁,
    idxTransport P F (idxPowLift' (F := F) U Y (e := k₁) (t := k₂ * k₁)
      (mul_comm k₂ k₁)) ψ₂, ?_, ?_⟩
  · exact (HomPf.mk_pow_lift' (F := F) rfl ψ₁).trans h₁
  · exact (HomPf.mk_pow_lift' (F := F) (mul_comm k₂ k₁) ψ₂).trans h₂

/-- ★★★★locator —— `Proposition 5.5, (ii)` の逆向きの辞書。 -/
def HomPf.mk_idxPow'.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 押し出しの添字での代表元と rootIso",
    sectionId := "frdi-prop-5-5" }

/-! ## ★10. 根がずれた合成の計算則

★★`Proposition 5.5, (ii)` の逆向きでは、根 `n` の対象から根 1 の対象への射に
`𝒞` の射の像を後合成する形が出る。
★`compRoot_eq_lift`(在庫)に `c = 1`、`PA = 1`、`PB = PE = n` を入れると、
**3 つの `rtRootIso` のうち 2 つが恒等になり**、残る 1 つだけが残る。 -/

/-- ★★★**`𝒞` の射の像を後合成する計算則**(根 `n` → 根 1 → 根 1)。 -/
theorem compRoot_toRootHom {A B E : C} {n : ℕ+}
    (f : HomRoot P F (⟨A, n⟩ : PfRootObj P F) ⟨B, 1⟩) (g : B ⟶ E) :
    compRoot P F f (toRootHom (F := F) g)
      = compPf P F f ((rtRootIso P F B E (show n = n * 1 from (mul_one n).symm)
          (show n = n * 1 from (mul_one n).symm)).inv (toRootHom (F := F) g)) := by
  have h := compRoot_eq_lift (F := F) f (toRootHom (F := F) g)
    (c := 1) (PA := 1) (PB := n) (PE := n)
    (hcA := by rw [one_mul, one_mul])
    (hcB := by rw [one_mul, one_mul])
    (hcE := by rw [one_mul, one_mul])
    (ef := 1) (eg := n) (er := 1)
    (hfA := (one_mul 1).symm) (hfB := (one_mul n).symm)
    (hgA := (mul_one n).symm) (hgE := (mul_one n).symm)
    (hrA := (one_mul 1).symm) (hrE := (one_mul n).symm)
  rw [h, rtRootIso_inv_eq_self, rtRootIso_hom_eq_self]

/-! ## ★11. 根の不変性と `toHomPf`

★★`rootIso` は**`frobTransport` を消す** —— これで
`rtRootIso` の `.inv` を `𝒞` の射の言葉に落とせる。 -/

/-- ★添字 `idxOne A B ⟶ (pushIdx a b).obj (idxOne A′ B′)`。 -/
noncomputable def idxOneHomPush {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) :
    idxOne P F A B ⟶ (pushIdx (F := F) a ha b hb hd).obj (idxOne P F A' B') :=
  Under.homMk
    (show (⟨(A, B)⟩ : BiFr P F) ⟶ (⟨(A', B')⟩ : BiFr P F) from ⟨(a, b), ha, hb, hd⟩)
    (WideSubcategory.hom_ext _
      (Prod.ext
        (show 𝟙 A ≫ a = a ≫ 𝟙 A' by rw [Category.id_comp, Category.comp_id])
        (show 𝟙 B ≫ b = b ≫ 𝟙 B' by rw [Category.id_comp, Category.comp_id])))

/-- ★★★**根の不変性の `toHomPf` での計算則** —— `.hom` は `frobTransport` を消す。 -/
theorem rootIso_hom_toHomPf {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) (φ : A ⟶ B) :
    (rootIso (F := F) a ha b hb hd).hom
        (toHomPf (F := F) (frobTransport (F := F) a ha b hb hd φ))
      = toHomPf (F := F) φ := by
  rw [toHomPf, rootIso_hom_mk]
  exact HomPf.mk_map (idxOneHomPush (F := F) a ha b hb hd) φ

/-- ★★★**`.inv` の計算則**。 -/
theorem rootIso_inv_toHomPf {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) (φ : A ⟶ B) :
    (rootIso (F := F) a ha b hb hd).inv (toHomPf (F := F) φ)
      = toHomPf (F := F) (frobTransport (F := F) a ha b hb hd φ) := by
  refine ((congrArg (rootIso (F := F) a ha b hb hd).inv
    (rootIso_hom_toHomPf a ha b hb hd φ)).symm).trans ?_
  exact Iso.hom_inv_id_apply (rootIso (F := F) a ha b hb hd) _

/-- ★★★★locator —— `Proposition 5.5, (ii)` の逆向きに要る合成の計算則。 -/
def compRoot_toRootHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 根がずれた合成の計算則",
    sectionId := "frdi-prop-5-5" }

/-! ## ★12. `idxPow` と押し出しは同じ元を与える

★★第 9 節で「2 つの添字対象は定義的に等しくない」と書いたが、
**余極限の元としては等しい**。★理由は 2 つ:
* 添字圏が**細い**(`idx_hom_ext`) —— `idxOne` からの 2 本の合成が一致する。
* `𝒞` が **totally epimorphic**(`rtExt` が epi) —— そこから脚が一致する。

★これで `compPf_mk_pair`(在庫、`idxMk` の形を要求する)が使えるようになる。 -/

/-- ★★★**`idxPow` と `idxPow'` は同じ元を与える**。 -/
theorem HomPf.mk_idxPow_eq (U V : C) (k : ℕ+) (ψ : rtObj P F U k ⟶ rtObj P F V k) :
    HomPf.mk (idxPow (F := F) U V k) ψ = HomPf.mk (idxPow' (F := F) U V k) ψ := by
  have hdd : P.degFr (rtExt P F U k) = P.degFr (rtExt P F V k) := by
    rw [rtExt_degFr, rtExt_degFr]
  set W := IsFiltered.max (idxPow (F := F) U V k) (idxPow' (F := F) U V k) with hW
  set u := IsFiltered.leftToMax (idxPow (F := F) U V k) (idxPow' (F := F) U V k) with hu
  set u' := IsFiltered.rightToMax (idxPow (F := F) U V k) (idxPow' (F := F) U V k) with hu'
  have hpar : (idxOneHom (F := F) (rtExt P F U k) (rtExt P F V k)
        (rtExt_frobType P F U k) (rtExt_frobType P F V k) hdd) ≫ u
      = (idxOneHomPush (F := F) (rtExt P F U k) (rtExt_frobType P F U k)
        (rtExt P F V k) (rtExt_frobType P F V k) hdd) ≫ u' :=
    idx_hom_ext _ _
  have h1 : rtExt P F U k ≫ u.right.hom.1 = rtExt P F U k ≫ u'.right.hom.1 :=
    congrArg (fun t : idxOne P F U V ⟶ W => t.right.hom.1) hpar
  have h2 : rtExt P F V k ≫ u.right.hom.2 = rtExt P F V k ≫ u'.right.hom.2 :=
    congrArg (fun t : idxOne P F U V ⟶ W => t.right.hom.2) hpar
  have hq1 : u.right.hom.1 = u'.right.hom.1 :=
    haveI := P.totEpiC U (rtObj P F U k) (rtExt P F U k)
    (cancel_epi (rtExt P F U k)).mp h1
  have hq2 : u.right.hom.2 = u'.right.hom.2 :=
    haveI := P.totEpiC V (rtObj P F V k) (rtExt P F V k)
    (cancel_epi (rtExt P F V k)).mp h2
  refine HomColim.sound _ W u u' (congrArg ULift.up ?_)
  refine frobTransport_eq (F := F) _ _ _ _ _ ψ _ ?_
  exact ((congrArg (fun t => ψ ≫ t) hq2).trans
    (frobTransport_spec (F := F) u'.right.hom.1 u'.right.property.1
      u'.right.hom.2 u'.right.property.2.1 u'.right.property.2.2 ψ)).trans
    (congrArg (fun t => t ≫ idxTransport P F u' ψ) hq1.symm)

/-- ★★★**辞書(`idxMk` の形)** —— `idxPow` の代表元も `rootIso` で `toHomPf` に戻る。 -/
theorem HomPf.mk_idxPow_rootIso (U V : C) (k : ℕ+) (ψ : rtObj P F U k ⟶ rtObj P F V k) :
    HomPf.mk (idxPow (F := F) U V k) ψ
      = (rootIso (F := F) (rtExt P F U k) (rtExt_frobType P F U k)
          (rtExt P F V k) (rtExt_frobType P F V k)
          (by rw [rtExt_degFr, rtExt_degFr])).hom (toHomPf (F := F) ψ) :=
  (HomPf.mk_idxPow_eq U V k ψ).trans (HomPf.mk_idxPow' U V k ψ)

/-- ★★★locator —— `Proposition 5.5, (ii)` の添字の同一視。 -/
def HomPf.mk_idxPow_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — idxPow と押し出しは同じ元を与える",
    sectionId := "frdi-prop-5-5" }

/-! ## ★13. 代表元に `𝒞` の射の像を後合成する

★★`toHomPf θ` を添字 `idxMk b e` へ押し上げる(`idxOneHom`)と、
`compPf_mk_pair`(在庫)がそのまま当たる。 -/

/-- ★★**`toHomPf` を添字 `idxMk a e` へ押し上げる**。 -/
theorem toHomPf_eq_mk_idxMk {A E A' E' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (e : E ⟶ E') (he : IsFrobeniusType P e) (hae : P.degFr a = P.degFr e) (ξ : A ⟶ E) :
    toHomPf (F := F) ξ
      = HomPf.mk (idxMk (P := P) (F := F) a e ha he hae)
        (frobTransport (F := F) a ha e he hae ξ) :=
  (HomPf.mk_map (idxOneHom (F := F) a e ha he hae) ξ).symm

/-- ★★★**代表元に `𝒞` の射の像を後合成する計算則**。 -/
theorem compPf_mk_toHomPf {A B E A' B' E' : C} (a : A ⟶ A') (b : B ⟶ B')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (φ : A' ⟶ B') (θ : B ⟶ E) (e : E ⟶ E') (he : IsFrobeniusType P e)
    (hbe : P.degFr b = P.degFr e) :
    compPf P F (HomPf.mk (idxMk (P := P) (F := F) a b ha hb hd) φ) (toHomPf (F := F) θ)
      = HomPf.mk (idxMk (P := P) (F := F) a e ha he (hd.trans hbe))
        (φ ≫ frobTransport (F := F) b hb e he hbe θ) := by
  rw [toHomPf_eq_mk_idxMk b hb e he hbe θ, compPf_mk_pair]
  rfl

/-- ★★★locator —— `Proposition 5.5, (ii)` の後合成の計算則。 -/
def compPf_mk_toHomPf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 代表元への後合成の計算則",
    sectionId := "frdi-prop-5-5" }

/-! ## ★14. 残る図式追跡の材料

★★第 21 弾で書き下した等式の材料をここに置く。
★要点は 3 つ:
* `kappaRep A n = rtLift A (n = n·1)` —— **`κ` の代表元は `rtLift` そのもの**。
* `θ_A = rtExt X k ≫ rtExt E n` —— **`θ_A` は `rtExt` の合成そのもの**。
* `frobTransport α β (α ≫ ξ) = ξ ≫ β` —— **前合成した射の遷移は後合成**。 -/

/-- ★★**`κ` の代表元は `rtLift` そのもの**。 -/
theorem kappaRep_eq_rtLift (A : C) (n : ℕ+) :
    kappaRep (P := P) (F := F) A n = rtLift P F A (show n = n * 1 from (mul_one n).symm) := by
  refine rtLift_uniq P F _ _ ?_ (rtLift_ext P F A _)
  show rtExt P F A 1 ≫ (rtOneInv (P := P) (F := F) A ≫ rtExt P F A n) = rtExt P F A n
  rw [← Category.assoc]
  have h : rtExt P F A 1 ≫ rtOneInv (P := P) (F := F) A = 𝟙 A :=
    @IsIso.hom_inv_id _ _ _ _ (rtExt P F A 1) (isIso_rtExt_one P F A)
  rw [h, Category.id_comp]

/-- ★★**`α` を前合成した射の遷移は「後合成」になる**。 -/
theorem frobTransport_comp_left {A' A'' B' B'' : C} (α : A' ⟶ A'') (hα : IsFrobeniusType P α)
    (β : B' ⟶ B'') (hβ : IsFrobeniusType P β) (hd : P.degFr α = P.degFr β)
    (ξ : A'' ⟶ B') :
    frobTransport (F := F) α hα β hβ hd (α ≫ ξ) = ξ ≫ β :=
  frobTransport_eq (F := F) α hα β hβ hd (α ≫ ξ) (ξ ≫ β) (Category.assoc _ _ _)

/-- ★★★**`θ_A` の正体** —— `rtExt` の合成そのもの。 -/
theorem theta_A_eq (A : C) (n k : ℕ+) :
    frobTransport (F := F)
        (rtLift P F A (show n = n * 1 from (mul_one n).symm))
        (rtLift_frobType P F A _)
        (rtLift P F (rtObj P F (rtObj P F A n) k) (show n = n * 1 from (mul_one n).symm))
        (rtLift_frobType P F _ _)
        (by rw [rtLift_degFr, rtLift_degFr])
        (rtOneInv (P := P) (F := F) A
          ≫ (rtExt P F A n ≫ rtExt P F (rtObj P F A n) k)
          ≫ rtExt P F (rtObj P F (rtObj P F A n) k) 1)
      = rtExt P F (rtObj P F A n) k ≫ rtExt P F (rtObj P F (rtObj P F A n) k) n := by
  refine frobTransport_eq (F := F) _ _ _ _ _ _ _ ?_
  have hk : rtOneInv (P := P) (F := F) A ≫ rtExt P F A n
      = rtLift P F A (show n = n * 1 from (mul_one n).symm) := kappaRep_eq_rtLift A n
  have hE : rtExt P F (rtObj P F (rtObj P F A n) k) 1
        ≫ rtLift P F (rtObj P F (rtObj P F A n) k) (show n = n * 1 from (mul_one n).symm)
      = rtExt P F (rtObj P F (rtObj P F A n) k) n := rtLift_ext P F _ _
  calc (rtOneInv (P := P) (F := F) A ≫ (rtExt P F A n ≫ rtExt P F (rtObj P F A n) k)
          ≫ rtExt P F (rtObj P F (rtObj P F A n) k) 1)
        ≫ rtLift P F (rtObj P F (rtObj P F A n) k) (show n = n * 1 from (mul_one n).symm)
      = (rtOneInv (P := P) (F := F) A ≫ rtExt P F A n) ≫ rtExt P F (rtObj P F A n) k
          ≫ (rtExt P F (rtObj P F (rtObj P F A n) k) 1
            ≫ rtLift P F (rtObj P F (rtObj P F A n) k)
              (show n = n * 1 from (mul_one n).symm)) := by simp only [Category.assoc]
    _ = rtLift P F A (show n = n * 1 from (mul_one n).symm)
        ≫ rtExt P F (rtObj P F A n) k ≫ rtExt P F (rtObj P F (rtObj P F A n) k) n := by
        rw [hk, hE]

/-- ★★★`θ_A` を遷移した先 —— `rtExt E n ≫ e`。 -/
theorem frobTransport_theta_A (A : C) (n k : ℕ+) {E' : C}
    (e : rtObj P F (rtObj P F (rtObj P F A n) k) n ⟶ E') (he : IsFrobeniusType P e)
    (hd : P.degFr (rtExt P F (rtObj P F A n) k) = P.degFr e) :
    frobTransport (F := F) (rtExt P F (rtObj P F A n) k)
        (rtExt_frobType P F _ _) e he hd
        (rtExt P F (rtObj P F A n) k ≫ rtExt P F (rtObj P F (rtObj P F A n) k) n)
      = rtExt P F (rtObj P F (rtObj P F A n) k) n ≫ e :=
  frobTransport_comp_left _ _ _ he hd _

/-- ★★★locator —— `Proposition 5.5, (ii)` の図式追跡の材料。 -/
def theta_A_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 辞書の図式追跡の材料",
    sectionId := "frdi-prop-5-5" }

/-! ## ★15. ★★★★★★★辞書(i′) —— 逆向きの核心

★★第 21 弾で書き下した等式を、第 14 節の材料で閉じる。
★中身は 5 段:
1. `compRoot_toRootHom` で `compRoot` を `compPf` に落とす（`rtRootIso` が 2 つ消える）
2. `rtRootIso_inv_toHomPf` で残る `rtRootIso` の `.inv` を `frobTransport` に落とす
3. `theta_A_eq` で左の `frobTransport` を `rtExt` の合成に潰す
4. `kappaRep_eq_rtLift` ＋ `frobTransport_spec` で右を潰す
5. `compPf_mk_final`（後合成の計算則 2 回）で両辺を突き合わせる -/

theorem inv_rtExt_one_eq (B : C) :
    @inv _ _ _ _ (rtExt P F B 1) (isIso_rtExt_one P F B) = rtOneInv (P := P) (F := F) B := rfl

theorem rtRootIso_inv_toHomPf (A B : C) {dA dB e tA tB : ℕ+} (hA : tA = e * dA)
    (hB : tB = e * dB) (φ : rtObj P F A dA ⟶ rtObj P F B dB) :
    (rtRootIso P F A B hA hB).inv (toHomPf (F := F) φ)
      = toHomPf (F := F) (frobTransport (F := F) (rtLift P F A hA) (rtLift_frobType P F A hA)
          (rtLift P F B hB) (rtLift_frobType P F B hB)
          (by rw [rtLift_degFr, rtLift_degFr]) φ) :=
  rootIso_inv_toHomPf _ _ _ _ _ φ

set_option maxHeartbeats 800000 in
/-- ★★★★**両辺を同じ添字に揃えて突き合わせる段**。 -/
theorem compPf_mk_final (A'' A : C) (n k : ℕ+)
    (ψ₁ : rtObj P F (rtObj P F A'' 1) k ⟶ rtObj P F (rtObj P F A n) k)
    (hdeg1 : P.degFr (rtExt P F (rtObj P F A'' 1) k)
      = P.degFr (rtExt P F (rtObj P F A n) k)) :
    compPf P F (HomPf.mk (idxMk (P := P) (F := F) (rtExt P F (rtObj P F A'' 1) k)
          (rtExt P F (rtObj P F A n) k) (rtExt_frobType P F _ _) (rtExt_frobType P F _ _)
          hdeg1) ψ₁)
        (toHomPf (F := F) (rtExt P F (rtObj P F A n) k
          ≫ rtExt P F (rtObj P F (rtObj P F A n) k) n))
      = toHomPf (F := F) (rtExt P F (rtObj P F A'' 1) k ≫ ψ₁
          ≫ rtExt P F (rtObj P F (rtObj P F A n) k) n) := by
  have hdeg2 : P.degFr (rtExt P F (rtObj P F A n) k)
      = P.degFr (rtExt P F (rtObj P F (rtObj P F (rtObj P F A n) k) n) k) := by
    rw [rtExt_degFr, rtExt_degFr]
  have hdeg3 : P.degFr (rtExt P F (rtObj P F A'' 1) k)
      = P.degFr (rtExt P F (rtObj P F (rtObj P F (rtObj P F A n) k) n) k) := hdeg1.trans hdeg2
  have hpush : toHomPf (F := F) (rtExt P F (rtObj P F A'' 1) k ≫ ψ₁
        ≫ rtExt P F (rtObj P F (rtObj P F A n) k) n)
      = HomPf.mk (idxMk (P := P) (F := F) (rtExt P F (rtObj P F A'' 1) k)
          (rtExt P F (rtObj P F (rtObj P F (rtObj P F A n) k) n) k)
          (rtExt_frobType P F _ _) (rtExt_frobType P F _ _) hdeg3)
        (frobTransport (F := F) (rtExt P F (rtObj P F A'' 1) k) (rtExt_frobType P F _ _)
          (rtExt P F (rtObj P F (rtObj P F (rtObj P F A n) k) n) k) (rtExt_frobType P F _ _)
          hdeg3 (rtExt P F (rtObj P F A'' 1) k ≫ ψ₁
            ≫ rtExt P F (rtObj P F (rtObj P F A n) k) n)) :=
    (HomPf.mk_map (idxOneHom (F := F) _ _ (rtExt_frobType P F _ _)
      (rtExt_frobType P F _ _) hdeg3) _).symm
  have hfin : frobTransport (F := F) (rtExt P F (rtObj P F A'' 1) k)
        (rtExt_frobType P F _ _)
        (rtExt P F (rtObj P F (rtObj P F (rtObj P F A n) k) n) k) (rtExt_frobType P F _ _)
        hdeg3 (rtExt P F (rtObj P F A'' 1) k ≫ (ψ₁
          ≫ rtExt P F (rtObj P F (rtObj P F A n) k) n))
      = (ψ₁ ≫ rtExt P F (rtObj P F (rtObj P F A n) k) n)
        ≫ rtExt P F (rtObj P F (rtObj P F (rtObj P F A n) k) n) k :=
    frobTransport_comp_left _ _ _ _ _ _
  rw [compPf_mk_toHomPf (P := P) (F := F) _ _ _ _ hdeg1 ψ₁ _
    (rtExt P F (rtObj P F (rtObj P F (rtObj P F A n) k) n) k)
    (rtExt_frobType P F _ _) hdeg2]
  rw [frobTransport_comp_left, hpush]
  exact congrArg (HomPf.mk _) ((Category.assoc _ _ _).symm.trans hfin.symm)

set_option maxHeartbeats 800000 in
/-- ★★★★★★★**辞書(i′)** —— `Proposition 5.5, (ii)` の逆向きの核心。

`idxPow` の添字での代表元に `𝒞` の Frobenius 拡大の像を後合成したものは、
**`κ` を前合成した形**に書き換わる。 -/
theorem compRoot_dictionary (A'' A : C) (n k : ℕ+)
    (ψ₁ : rtObj P F (rtObj P F A'' 1) k ⟶ rtObj P F (rtObj P F A n) k) :
    compRoot P F
        (show HomRoot P F (⟨A'', n⟩ : PfRootObj P F) ⟨A, 1⟩ from
          HomPf.mk (idxPow (F := F) (rtObj P F A'' 1) (rtObj P F A n) k) ψ₁)
        (toRootHom (F := F) (rtExt P F A n ≫ rtExt P F (rtObj P F A n) k))
      = compRoot P F (pfKappa (F := F) A'' n)
        (toRootHom (F := F)
          ((rtExt P F A'' 1 ≫ rtExt P F (rtObj P F A'' 1) k) ≫ ψ₁)) := by
  have hdeg1 : P.degFr (rtExt P F (rtObj P F A'' 1) k)
      = P.degFr (rtExt P F (rtObj P F A n) k) := by rw [rtExt_degFr, rtExt_degFr]
  have hAAinv : rtOneInv (P := P) (F := F) A'' ≫ rtExt P F A'' 1 = 𝟙 (rtObj P F A'' 1) :=
    @IsIso.inv_hom_id _ _ _ _ (rtExt P F A'' 1) (isIso_rtExt_one P F A'')
  have hinput : rtOneInv (P := P) (F := F) A''
        ≫ ((rtExt P F A'' 1 ≫ rtExt P F (rtObj P F A'' 1) k) ≫ ψ₁)
        ≫ rtExt P F (rtObj P F (rtObj P F A n) k) 1
      = rtExt P F (rtObj P F A'' 1) k ≫ ψ₁
        ≫ rtExt P F (rtObj P F (rtObj P F A n) k) 1 := by
    simp only [← Category.assoc, hAAinv]
    simp
  have hkappa : kappaRep (P := P) (F := F) A'' n
        ≫ frobTransport (F := F)
          (rtLift P F A'' (show n = n * 1 from (mul_one n).symm))
          (rtLift_frobType P F A'' _)
          (rtLift P F (rtObj P F (rtObj P F A n) k)
            (show n = n * 1 from (mul_one n).symm))
          (rtLift_frobType P F _ _) (by rw [rtLift_degFr, rtLift_degFr])
          (rtExt P F (rtObj P F A'' 1) k ≫ ψ₁
            ≫ rtExt P F (rtObj P F (rtObj P F A n) k) 1)
      = rtExt P F (rtObj P F A'' 1) k ≫ ψ₁
        ≫ rtExt P F (rtObj P F (rtObj P F A n) k) n := by
    rw [kappaRep_eq_rtLift]
    rw [← frobTransport_spec (F := F)
      (rtLift P F A'' (show n = n * 1 from (mul_one n).symm)) (rtLift_frobType P F A'' _)
      (rtLift P F (rtObj P F (rtObj P F A n) k) (show n = n * 1 from (mul_one n).symm))
      (rtLift_frobType P F _ _) (by rw [rtLift_degFr, rtLift_degFr]) _]
    rw [Category.assoc, Category.assoc, rtLift_ext]
  rw [compRoot_toRootHom, compRoot_toRootHom, toRootHom, toRootHom,
    rtRootIso_inv_toHomPf, rtRootIso_inv_toHomPf, inv_rtExt_one_eq, inv_rtExt_one_eq,
    theta_A_eq, hinput, show pfKappa (F := F) A'' n
      = toHomPf (F := F) (kappaRep (P := P) (F := F) A'' n) from rfl,
    ← toHomPf_comp, hkappa]
  exact compPf_mk_final A'' A n k ψ₁ hdeg1

/-- ★★★★★★locator —— `Proposition 5.5, (ii)` の逆向きの辞書。 -/
def compRoot_dictionary.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 逆向きの辞書(κ を前合成した形へ)",
    sectionId := "frdi-prop-5-5" }

/-! ## ★16. `κ` に沿った持ち上げ(一般形)

★★`rootLift`(第 3 節)は始域が `⟨A,1⟩` の場合だったが、
**任意の始域 `X₀` から `⟨V,1⟩` への Frobenius 型射**についても同じことが成り立つ。
★存在は `pfRoot_frob_div`、同型性は `κ` と `g` がどちらも mono であることから。

★★★これが `Proposition 5.5, (ii)` の**全射性の三角形**を作る道具である ——
「`κ` の合成則」を別に用意しなくても、**持ち上げを方程式で直接定義すれば済む**。 -/

/-- ★★★★**`κ` に沿った持ち上げ(一般形)** —— 任意の始域から。 -/
theorem exists_kappaLift (hfi : IsOfFrobeniusIsotropicType P) {X₀ : PfRootObj P F} {V : C}
    (K : ℕ+) (g : X₀ ⟶ (⟨V, 1⟩ : PfRootObj P F))
    (hdeg : (pfRootPre P F).degFr g = K) :
    ∃ e : X₀ ⟶ (⟨V, K⟩ : PfRootObj P F), e ≫ pfKappa (F := F) V K = g := by
  refine pfRoot_frob_div (F := F) (n := K) (m := 1) hfi g (pfKappa (F := F) V K) ?_
    (pfKappa_frobType hfi V K) (pfKappa_degFr V K)
  rw [hdeg, mul_one]

/-- ★★★★**持ち上げは同型**(`g` が Frobenius 型なら)。 -/
theorem isIso_kappaLift (hfi : IsOfFrobeniusIsotropicType P) {X₀ : PfRootObj P F} {V : C}
    (K : ℕ+) (g : X₀ ⟶ (⟨V, 1⟩ : PfRootObj P F))
    (hg : IsFrobeniusType (pfRootPre P F) g) (hdeg : (pfRootPre P F).degFr g = K)
    (e : X₀ ⟶ (⟨V, K⟩ : PfRootObj P F)) (he : e ≫ pfKappa (F := F) V K = g) : IsIso e := by
  obtain ⟨e', he'⟩ := pfRoot_frob_div (F := F) (n := K) (m := 1) hfi
    (pfKappa (F := F) V K) g (by rw [pfKappa_degFr, mul_one]) hg hdeg
  haveI hm1 : Mono g := pfRoot_frobTypeMono hfi _ hg
  haveI hm2 : Mono (pfKappa (F := F) V K) := pfKappa_mono hfi V K
  refine ⟨e', ?_, ?_⟩
  · refine (cancel_mono g).mp ?_
    rw [Category.assoc, he', he, Category.id_comp]
  · refine (cancel_mono (pfKappa (F := F) V K)).mp ?_
    rw [Category.assoc, he, he', Category.id_comp]

/-- ★★★**`κ` に沿った標準の持ち上げ**。 -/
noncomputable def kappaLift (hfi : IsOfFrobeniusIsotropicType P) {X₀ : PfRootObj P F} {V : C}
    (K : ℕ+) (g : X₀ ⟶ (⟨V, 1⟩ : PfRootObj P F))
    (hdeg : (pfRootPre P F).degFr g = K) : X₀ ⟶ (⟨V, K⟩ : PfRootObj P F) :=
  (exists_kappaLift hfi K g hdeg).choose

@[simp] theorem kappaLift_spec (hfi : IsOfFrobeniusIsotropicType P) {X₀ : PfRootObj P F}
    {V : C} (K : ℕ+) (g : X₀ ⟶ (⟨V, 1⟩ : PfRootObj P F))
    (hdeg : (pfRootPre P F).degFr g = K) :
    kappaLift (F := F) hfi K g hdeg ≫ pfKappa (F := F) V K = g :=
  (exists_kappaLift hfi K g hdeg).choose_spec

theorem kappaLift_isIso (hfi : IsOfFrobeniusIsotropicType P) {X₀ : PfRootObj P F} {V : C}
    (K : ℕ+) (g : X₀ ⟶ (⟨V, 1⟩ : PfRootObj P F))
    (hg : IsFrobeniusType (pfRootPre P F) g) (hdeg : (pfRootPre P F).degFr g = K) :
    IsIso (kappaLift (F := F) hfi K g hdeg) :=
  isIso_kappaLift hfi K g hg hdeg _ (kappaLift_spec hfi K g hdeg)

/-- ★★**`κ` と合成して等しければ射そのものが等しい**。 -/
theorem kappaLift_ext (hfi : IsOfFrobeniusIsotropicType P) {X₀ : PfRootObj P F} {V : C}
    {K : ℕ+} (e e' : X₀ ⟶ (⟨V, K⟩ : PfRootObj P F))
    (h : e ≫ pfKappa (F := F) V K = e' ≫ pfKappa (F := F) V K) : e = e' := by
  haveI := pfKappa_mono (F := F) hfi V K
  exact (cancel_mono (pfKappa (F := F) V K)).mp h

/-- ★★★★locator —— `Proposition 5.5, (ii)` の全射性に要る持ち上げ。 -/
def kappaLift.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — κ に沿った持ち上げ(一般形)",
    sectionId := "frdi-prop-5-5" }

/-! ## ★17. ★★★★★★全射性の三角形

★★第 22 弾の辞書と第 23 弾の持ち上げを合わせると、
`Proposition 5.5, (ii)` の全射性に要る三角形が**ただちに出る**。

★★同じ 1 本が `A` の側（構造射）にも `B` の側（値）にも当たる ——
`e_β` は `A''`・`n`・`k` にしか依らないからである。 -/

set_option maxHeartbeats 800000 in
/-- ★★★★★★**全射性の三角形** —— 辞書からただちに出る。

`ε ≫ iA_W = e_β ≫ rootMap ψ₁ K`（`K = k·n`）。
★★`κ` が mono なので、`κ` と合成した式を比べればよい:
* 左 —— `rootLift_spec` で `[rtExt A n ≫ rtExt X k]` になり、**辞書**で
  `pfKappa A'' n ≫ [β ≫ ψ₁]` になる。
* 右 —— `rootMap_spec` と `kappaLift_spec` で `(pfKappa A'' n ≫ [β]) ≫ [ψ₁]` になる。 -/
theorem surj_triangle (hfi : IsOfFrobeniusIsotropicType P) (A'' A : C) (n k : ℕ+)
    (ψ₁ : rtObj P F (rtObj P F A'' 1) k ⟶ rtObj P F (rtObj P F A n) k)
    (hA : P.degFr (rtExt P F A n ≫ rtExt P F (rtObj P F A n) k) = k * n)
    (hB : (pfRootPre P F).degFr (pfKappa (F := F) A'' n
      ≫ toRootHom (F := F) (rtExt P F A'' 1 ≫ rtExt P F (rtObj P F A'' 1) k)) = k * n) :
    (show HomRoot P F (⟨A'', n⟩ : PfRootObj P F) ⟨A, 1⟩ from
        HomPf.mk (idxPow (F := F) (rtObj P F A'' 1) (rtObj P F A n) k) ψ₁)
        ≫ rootLift (F := F) hfi
          (rtExt P F A n ≫ rtExt P F (rtObj P F A n) k) (k * n) hA
      = kappaLift (F := F) hfi (k * n)
          (pfKappa (F := F) A'' n
            ≫ toRootHom (F := F) (rtExt P F A'' 1 ≫ rtExt P F (rtObj P F A'' 1) k)) hB
        ≫ rootMap (F := F) hfi ψ₁ (k * n) := by
  refine kappaLift_ext hfi _ _ ?_
  rw [Category.assoc, rootLift_spec, Category.assoc, rootMap_spec, ← Category.assoc,
    kappaLift_spec, Category.assoc]
  exact (compRoot_dictionary A'' A n k ψ₁).trans
    (congrArg (fun t : (⟨A'', 1⟩ : PfRootObj P F) ⟶ ⟨rtObj P F (rtObj P F A n) k, 1⟩ =>
      compRoot P F (pfKappa (F := F) A'' n) t) (toRootHom_comp _ _))

/-- ★★★★★★locator —— `Proposition 5.5, (ii)` の全射性の三角形。 -/
def surj_triangle.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 全射性の三角形",
    sectionId := "frdi-prop-5-5" }

/-! ## ★18. 三角形の一般形と次数の計算 -/

set_option maxHeartbeats 800000 in
/-- ★★★★★★**全射性の三角形**(`K` を一般にした版)。

★`K` を仮引数にしておくと、`biratPfDeg W`(定義上は `P.degFr W₁`)を
そのまま入れられる —— 添字の根がずれずに済む。 -/
theorem surj_triangle' (hfi : IsOfFrobeniusIsotropicType P) (A'' A : C) (n k K : ℕ+)
    (ψ₁ : rtObj P F (rtObj P F A'' 1) k ⟶ rtObj P F (rtObj P F A n) k)
    (hA : P.degFr (rtExt P F A n ≫ rtExt P F (rtObj P F A n) k) = K)
    (hB : (pfRootPre P F).degFr (pfKappa (F := F) A'' n
      ≫ toRootHom (F := F) (rtExt P F A'' 1 ≫ rtExt P F (rtObj P F A'' 1) k)) = K) :
    (show HomRoot P F (⟨A'', n⟩ : PfRootObj P F) ⟨A, 1⟩ from
        HomPf.mk (idxPow (F := F) (rtObj P F A'' 1) (rtObj P F A n) k) ψ₁)
        ≫ rootLift (F := F) hfi
          (rtExt P F A n ≫ rtExt P F (rtObj P F A n) k) K hA
      = kappaLift (F := F) hfi K
          (pfKappa (F := F) A'' n
            ≫ toRootHom (F := F) (rtExt P F A'' 1 ≫ rtExt P F (rtObj P F A'' 1) k)) hB
        ≫ rootMap (F := F) hfi ψ₁ K := by
  refine kappaLift_ext hfi _ _ ?_
  rw [Category.assoc, rootLift_spec, Category.assoc, rootMap_spec, ← Category.assoc,
    kappaLift_spec, Category.assoc]
  exact (compRoot_dictionary A'' A n k ψ₁).trans
    (congrArg (fun t : (⟨A'', 1⟩ : PfRootObj P F) ⟶ ⟨rtObj P F (rtObj P F A n) k, 1⟩ =>
      compRoot P F (pfKappa (F := F) A'' n) t) (toRootHom_comp _ _))

/-- ★次数の計算 —— `K = k · n`(`A` の側)。 -/
theorem surj_degA (A : C) (n k : ℕ+) :
    P.degFr (rtExt P F A n ≫ rtExt P F (rtObj P F A n) k) = k * n := by
  rw [P.degFr_comp, rtExt_degFr, rtExt_degFr]

/-- ★次数の計算 —— `K = k · n`(`κ ≫ [β]` の側)。 -/
theorem surj_degB (A'' : C) (n k : ℕ+) :
    (pfRootPre P F).degFr (pfKappa (F := F) A'' n
      ≫ toRootHom (F := F) (rtExt P F A'' 1 ≫ rtExt P F (rtObj P F A'' 1) k)) = k * n := by
  rw [(pfRootPre P F).degFr_comp, toRootHom_degFr', pfKappa_degFr, P.degFr_comp,
    rtExt_degFr, rtExt_degFr, mul_one]

/-- ★★★★★★locator —— `Proposition 5.5, (ii)` の全射性の三角形(一般形)。 -/
def surj_triangle'.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 全射性の三角形(K を一般にした版)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
