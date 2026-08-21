/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm62Psi
import ABC3.Found.FrdI.Thm52Frob
import ABC3.Found.FrdI.Prop21
import ABC3.Found.FrdI.PfGp
import ABC3.Found.FrdI.ModelRigid

/-!
# [FrdI] Theorem 6.2, (ii) —— Frobenius が定める `Ψ` は naive Frobenius 関手

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.111。

原文 (FrdI p.111):
> which is isomorphic to the naive Frobenius functor [of degree p on C] of Propo-

## ★★原文が一言で畳んだもの

原文の証明は「in a word, that both functors are obtained by raising to the p-th power」
の一言である。★その中身は 2 つに割れる:

1. **一般形** —— `𝟭 𝒞 ⟶ G` の成分がすべて**次数 `d` の Frobenius 型**なら、
   `G ≅ Ψ_naive`(`Proposition 2.1, (i)`)。
   ★これは `Definition 1.3, (ii)` の**本質的一意性**(`frobDegUniq`)そのものであり、
   自然性は `Proposition 1.10, (i)` の**一意性**から出る ——
   在庫の `nfCompare` の議論を「もう 1 つの選択」から「任意の関手」へ広げただけである。
2. **標数 `p` の場合の当てはめ** —— Frobenius は因子を `p` 倍し、有理函数を `p` 乗する。
   ★model Frobenioid の言葉では `Φ` も `B` も **`p` 倍**であり(`frobModelHom`)、
   `A ⟶ Ψ(A)` は `(𝟙, 0, p, 0)`(`frobUnitApp`)である。
   ★★**それが Frobenius 型であることは `model_frobeniusType_iff`(在庫)で
   「`Div = 0` かつ底が同型」に落ちる** —— model Frobenioid では
   **すべての射が co-angular**(`model_coAngular`)だからである。

## ★本ファイルで閉じること

| 宣言 | 中身 |
|---|---|
| `frobCompare` / `_isIso` / `_sq` / `_naturality` | ★一般形の比較射 |
| `naiveFrobIso_of_unit` | ★★**`G ≅ Ψ_naive`** |
| `gpMap_nsmulHom` | `gpMap` は `n` 倍と交換する |
| `frobModelHom` | ★Frobenius が定める `ModelData` の射(`Φ`・`B` を `p` 倍) |
| `frobUnitApp` | `A ⟶ Ψ(A)` |
| `thm62FrobIso` | ★★★**[FrdI] Theorem 6.2, (ii)** |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

/-! ## ★1. 一般形 —— 次数 `d` の Frobenius 型の単位を持つ関手は naive Frobenius -/

section NaiveFrobUniq

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P) (d : ℕ+)
  {G : C ⥤ C} (a : ∀ A : C, A ⟶ G.obj A)
  (hft : ∀ A : C, IsFrobeniusType P (a A)) (hdeg : ∀ A : C, P.degFr (a A) = d)

/-- ★**比較射** —— `Definition 1.3, (ii)` の本質的一意性が与える同型。 -/
noncomputable def frobCompare (A : C) : G.obj A ⟶ nfObj P F d A :=
  (F.frobDegUniq A (G.obj A) (nfObj P F d A) (a A) (nfHom P F d A)
    (hft A) (nfHom_frobType P F d A) (by rw [hdeg A, nfHom_degFr])).choose

theorem frobCompare_isIso (A : C) : IsIso (frobCompare P F d a hft hdeg A) :=
  (F.frobDegUniq A (G.obj A) (nfObj P F d A) (a A) (nfHom P F d A)
    (hft A) (nfHom_frobType P F d A) (by rw [hdeg A, nfHom_degFr])).choose_spec.1

theorem frobCompare_sq (A : C) :
    a A ≫ frobCompare P F d a hft hdeg A = nfHom P F d A :=
  (F.frobDegUniq A (G.obj A) (nfObj P F d A) (a A) (nfHom P F d A)
    (hft A) (nfHom_frobType P F d A) (by rw [hdeg A, nfHom_degFr])).choose_spec.2

/-- ★★**比較射は自然** —— `Proposition 1.10, (i)` の**一意性**だけから出る。 -/
theorem frobCompare_naturality
    (hnat : ∀ {A B : C} (φ : A ⟶ B), φ ≫ a B = a A ≫ G.map φ) {A B : C} (φ : A ⟶ B) :
    G.map φ ≫ frobCompare P F d a hft hdeg B
      = frobCompare P F d a hft hdeg A ≫ nfMap P F d φ := by
  refine prop_1_10_i_uniq P.totEpiC φ (a A) (nfHom P F d B) _ _ ?_ ?_
  · rw [← Category.assoc, ← hnat φ, Category.assoc, frobCompare_sq]
  · rw [← Category.assoc, frobCompare_sq, ← nfMap_sq P F d φ]

/-- ★★★★**次数 `d` の Frobenius 型の単位を持つ関手は naive Frobenius 関手と同型**。

★`Proposition 2.1, (i)` の「well-defined up to isomorphism of functors」を
**任意の関手へ広げた形**である。 -/
noncomputable def naiveFrobIso_of_unit
    (hnat : ∀ {A B : C} (φ : A ⟶ B), φ ≫ a B = a A ≫ G.map φ) :
    G ≅ naiveFrob P F d :=
  NatIso.ofComponents
    (fun A => @asIso _ _ _ _ (frobCompare P F d a hft hdeg A)
      (frobCompare_isIso P F d a hft hdeg A))
    (fun φ => frobCompare_naturality P F d a hft hdeg hnat φ)

end NaiveFrobUniq

/-! ## ★2. 標数 `p` の Frobenius —— `Φ` も `B` も `p` 倍 -/

section FrobModel

variable {D : Type u} [Category.{v} D] {M : ModelData.{v, u, w} D}

/-- ★`gpMap` は `n` 倍と交換する(`gp_hom_ext` で生成元だけ見る)。 -/
theorem gpMap_nsmulHom (N : Type w) [AddCommMonoid N] (n : ℕ) :
    gpMap N (nsmulHom n N) = nsmulHom n (Gp N) :=
  gp_hom_ext (fun m => by
    rw [gpMap_toGp]
    show toGp N (n • m) = n • toGp N m
    exact map_nsmul (toGpHom N) n m)

/-- ★★**Frobenius が定める `ModelData` の射** —— 底は恒等、`Φ` も `B` も `p` 倍。

★これが原文の「both functors are obtained by raising to the p-th power」の
`ModelData` の側である。 -/
noncomputable def frobModelHom (M : ModelData.{v, u, w} D) (p : ℕ+) :
    ModelDataHomOver (𝟭 D) M M where
  phiHom d := nsmulHom ((p : ℕ+) : ℕ) (M.phi.val d)
  phiNat f x := (map_nsmul (M.phi.map f) _ x).symm
  bmonHom d := nsmulHom ((p : ℕ+) : ℕ) (M.bmon.val d)
  bmonNat f x := (map_nsmul (M.bmon.map f) _ x).symm
  divCompat d u := by
    show M.divB d (((p : ℕ+) : ℕ) • u) = gpMap _ (nsmulHom _ _) (M.divB d u)
    rw [gpMap_nsmulHom, map_nsmul]
    rfl

/-- ★★**単位 `A ⟶ Ψ(A)`** —— 底は恒等、次数 `p`、零因子も有理函数も `0`。 -/
noncomputable def frobUnitApp (p : ℕ+) (A : ModelData.Obj M) :
    A ⟶ (thm62Psi (frobModelHom M p)).obj A :=
  { base := 𝟙 A.base, div := 0, deg := p, u := 0
    cond := by
      show ((p : ℕ+) : ℕ) • A.cls + toGpHom _ (0 : M.phi.val A.base)
        = M.phi.gpMapOn (𝟙 A.base) (gpMap _ (nsmulHom ((p : ℕ+) : ℕ) _) A.cls)
          + M.divB _ (0 : M.bmon.val A.base)
      rw [map_zero, add_zero, map_zero, add_zero, ModelData.gpMapOn_id_apply,
        gpMap_nsmulHom]
      rfl }

/-- ★★**単位は Frobenius 型** —— model Frobenioid では
**すべての射が co-angular**(`model_coAngular`)なので、
`Div = 0` かつ底が同型を見れば足りる。 -/
theorem frobUnitApp_frobType (h : ModelData.Hyp M) (p : ℕ+) (A : ModelData.Obj M) :
    IsFrobeniusType (ModelData.modelPre h) (frobUnitApp p A) :=
  (ModelData.model_frobeniusType_iff h _).mpr ⟨rfl, by
    show IsIso (𝟙 A.base); infer_instance⟩

theorem frobUnitApp_degFr (h : ModelData.Hyp M) (p : ℕ+) (A : ModelData.Obj M) :
    (ModelData.modelPre h).degFr (frobUnitApp p A) = p := rfl

/-- ★★**単位は自然** —— 4 成分の計算。次数の成分だけ `mul_comm` が要る。 -/
theorem frobUnitApp_naturality (p : ℕ+) {A B : ModelData.Obj M} (φ : A ⟶ B) :
    φ ≫ frobUnitApp p B = frobUnitApp p A ≫ (thm62Psi (frobModelHom M p)).map φ := by
  refine ModelData.Hom.ext ?_ ?_ ?_ ?_
  · show φ.base ≫ 𝟙 B.base = 𝟙 A.base ≫ (𝟭 D).map φ.base
    rw [Category.comp_id, Category.id_comp]; rfl
  · show M.phi.map φ.base (0 : M.phi.val B.base) + ((p : ℕ+) : ℕ) • φ.div
      = M.phi.map (𝟙 A.base) (((p : ℕ+) : ℕ) • φ.div) + (φ.deg : ℕ) • (0 : M.phi.val A.base)
    rw [map_zero, zero_add, smul_zero, add_zero]
    exact (M.phi.map_id A.base _).symm
  · show p * φ.deg = φ.deg * p
    exact mul_comm _ _
  · show M.bmon.map φ.base (0 : M.bmon.val B.base) + ((p : ℕ+) : ℕ) • φ.u
      = M.bmon.map (𝟙 A.base) (((p : ℕ+) : ℕ) • φ.u) + (φ.deg : ℕ) • (0 : M.bmon.val A.base)
    rw [map_zero, zero_add, smul_zero, add_zero]
    exact (M.bmon.map_id A.base _).symm

/-- ★★★★★**[FrdI] Theorem 6.2, (ii)** —— Frobenius が定める `Ψ` は
`Proposition 2.1` の **naive Frobenius 関手**(次数 `p`)と同型。

原文 (FrdI p.111):
> which is isomorphic to the naive Frobenius functor [of degree p on C] of Propo- -/
noncomputable def thm62FrobIso (h : ModelData.Hyp M)
    (Fc : FrobenioidCore (ModelData.modelPre h)) (p : ℕ+) :
    thm62Psi (frobModelHom M p) ≅ naiveFrob (ModelData.modelPre h) Fc p :=
  naiveFrobIso_of_unit (ModelData.modelPre h) Fc p (frobUnitApp p)
    (frobUnitApp_frobType h p) (frobUnitApp_degFr h p)
    (fun φ => frobUnitApp_naturality p φ)

end FrobModel

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Theorem 6.2, (ii)`。 -/
def thm62FrobIso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (ii) — Frobenius が定める Ψ は naive Frobenius 関手と同型",
    sectionId := "frdi-thm-6-2" }

/-- ★★locator —— `Proposition 2.1, (i)` の一意性を任意の関手へ広げた形。 -/
def naiveFrobIso_of_unit.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 44,
    item := "Proposition 2.1, (i) — 次数 d の Frobenius 型の単位を持つ関手は Ψ_naive と同型",
    sectionId := "frdi-prop-2-1" }

end ABC3.Found.FrdI
