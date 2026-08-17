import ABC3.Found.FrdI.Prop33UnTr

/-!
# [FrdI] Proposition 3.3, (iv) の最終文 —— 射の類型が `𝒞^istr` から来ること(続き)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.60。

原文 (FrdI p.60):
> phism; morphism of a given Frobenius degree) if and only if it arises from such

原文 (FrdI p.60):
> an arrow of Cistr.

★原文は **9 クラス**を挙げる: Frobenius 型 / pre-step / base-isomorphism / 同型 /
pull-back / isometry / co-angular / LB-invertible / 与えられた Frobenius 次数。

## ★仕分けの続き(`UnTr.lean` で 5 クラス、ここで 3 クラス)

| 群 | クラス | 実装 |
|---|---|---|
| 不変量だけで決まる | 次数 / base-isomorphism / isometry / pre-step | `UnTr.lean`(`Iff.rfl`) |
| 圏論的だが取れる | 同型 | `UnTr.lean`(`unTr_isIso_iff`) |
| ★**両側で自明**になる | **co-angular** | ★**本ファイル** |
| ★**上の 2 群に潰れる** | **Frobenius 型** / **LB-invertible** | ★**本ファイル** |
| 残り | pull-back | ★`⟸` は `unTr_isPullBack_of_istr`(子)、`⟹` は未 |

★★**要点は 2 つ**:
1. ★`𝒞^istr` も `𝒞^un-tr` も**全対象が isotropic** なので、
   `Proposition 1.4, (i)` により**どちらでも全射が co-angular** ——
   co-angular の同値は**両側が真**という形で成り立つ
2. ★★`Proposition 1.4, (i)` の「In particular」により、isotropic 型では
   **Frobenius 型 = 等長 ∧ 底同型**に潰れる。★その 2 条件は
   `Iff.rfl` で商を渡る(`UnTr.lean`)ので、**Frobenius 型も渡る**
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

include P in
/-- ★★**`𝒞^istr` でもすべての射が co-angular** —— 全対象が isotropic だから。

★`𝒞^un-tr` 側は子の `unTr_coAngular`。★**両側が真**なので同値は自明に成り立つ。 -/
theorem istr_coAngular_all (Fc : FrobenioidCore P) {A B : Istr P} (f : A ⟶ B) :
    IsCoAngular (istrPre P Fc) f :=
  prop_1_4_i (istrPre P Fc) f (fun X _ => istr_isotropic P Fc X)

include P in
/-- ★★★**co-angular は商を渡る**(原文の 9 クラスのうち 1 つ)。

★★**両側が自明に真**という形で成り立つ —— `𝒞^istr` も `𝒞^un-tr` も
全対象が isotropic だからである。 -/
theorem unTr_isCoAngular_iff (Fc : FrobenioidCore P) {A B : Istr P} (α : A.obj ⟶ B.obj) :
    IsCoAngular (unTrPre P Fc)
        (toHomUnTr P α : (show UnTr P from A) ⟶ (show UnTr P from B))
      ↔ IsCoAngular (istrPre P Fc) (ObjectProperty.homMk α : A ⟶ B) :=
  ⟨fun _ => istr_coAngular_all P Fc _, fun _ => unTr_coAngular P Fc _⟩

include P in
/-- ★★★**Frobenius 型は商を渡る**。

★★`Proposition 1.4, (i)` の「In particular」で**両側とも
「等長 ∧ 底同型」に潰れる**(`unTr_isFrobeniusType_iff` / `prop_1_4_i_frobeniusType`)。
★その 2 条件は `Iff.rfl` で渡る(`unTr_isIsometric_iff` / `unTr_isBaseIsomorphism_iff`)。 -/
theorem unTr_isFrobeniusType_iff' (Fc : FrobenioidCore P) {A B : Istr P}
    (α : A.obj ⟶ B.obj) :
    IsFrobeniusType (unTrPre P Fc)
        (toHomUnTr P α : (show UnTr P from A) ⟶ (show UnTr P from B))
      ↔ IsFrobeniusType (istrPre P Fc) (ObjectProperty.homMk α : A ⟶ B) := by
  rw [unTr_isFrobeniusType_iff P Fc,
    prop_1_4_i_frobeniusType (istrPre P Fc) _ (fun X _ => istr_isotropic P Fc X)]
  exact Iff.rfl

include P in
/-- ★★★**LB-invertible は商を渡る**。

★`IsLBInvertible = IsCoAngular ∧ IsIsometric` で、
co-angular は両側で自明、isometry は `Iff.rfl`。 -/
theorem unTr_isLBInvertible_iff (Fc : FrobenioidCore P) {A B : Istr P}
    (α : A.obj ⟶ B.obj) :
    IsLBInvertible (unTrPre P Fc)
        (toHomUnTr P α : (show UnTr P from A) ⟶ (show UnTr P from B))
      ↔ IsLBInvertible (istrPre P Fc) (ObjectProperty.homMk α : A ⟶ B) :=
  ⟨fun h => ⟨istr_coAngular_all P Fc _, h.2⟩,
   fun h => ⟨unTr_coAngular P Fc _, h.2⟩⟩

/-! ## ★残る 1 クラス —— pull-back の `⟹`

★`⟸`(`𝒞^istr` の pull-back は商へ渡る)は子が取った(`unTr_isPullBack_of_istr`)。
★★**`⟹` は未実装**である —— `𝒞^un-tr` で pull-back なら
代表元が `𝒞^istr` で pull-back か、を言う必要がある。
★商の普遍性が `𝒞^istr` の普遍性を復元するか、が要点で、
**単位同値のぶんの自由度をどう潰すか**が残っている。
-/

end ABC3.Found.FrdI
