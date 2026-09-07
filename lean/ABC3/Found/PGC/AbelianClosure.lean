import ABC3.Found.PGC.AbelianClosureSplit
import ABC3.Found.PGC.ProfiniteUnitsTorsion
import ABC3.Found.PGC.QpNonAbelian
import Mathlib.FieldTheory.Galois.Abelian

/-!
# `K^ab`(最大アーベル拡大)と `K^LT ≤ K^ab`(`sorry` 無し)

[Yoshida08] Theorem 6.15(Local Kronecker-Weber)へ向かう決定 D29 の道 B の
節点 **B1**。原典が記号 `K^ab` を説明なしに使うので、その**語彙を用意する**のが
本ファイルの仕事である。

```
K^ab := (K̄)^{⁅Γ_K, Γ_K⁆‾}   (Γ_K = Gal(K̄/K)、‾ は位相的閉包)
```

## ★到達点

* `abelianClosure`——`K^ab` を中間体として定義する(★木にも mathlib にも無かった)
* `instIsGaloisAbelianClosure` / `instNormalAbelianClosure`——`K^ab/K` は Galois
* `instIsAbelianGaloisAbelianClosure`——★`Gal(K^ab/K)` は**可換**
  (mathlib の `IsAbelianGalois` で言う)
* `le_abelianClosure_iff`——★**`K^ab` は最大アーベル部分拡大である**
  (`K` 上正規な `A` について `A ≤ K^ab ↔ Gal(A/K)` が可換)
* `unramifiedClosure_le_abelianClosure`(`K^ur ≤ K^ab`)
* `lubinTateClosure_le_abelianClosure`(`K_π ≤ K^ab`)
* `lubinTateClosure_sup_unramifiedClosure_le_abelianClosure`
  および仮説なしの `exists_lubinTate_le_abelianClosure`——**`K^LT ≤ K^ab`**

★★**`≤` であって `=` ではない。** 逆向き `K^ab ≤ K^LT` は B5(`E_σ ⊆ K_π`)を
経由して B6 で初めて出る。本ファイルはその半分しか主張しない。

## ★設計——抽象核と具体層

`## 0` の 8 本は**分岐・付値・Lubin-Tate の語彙が 1 語も出ない**
(一般の Galois 拡大 `E/k` と、その中間体だけ)。具体層(`## 1` 以降)は
そこへ `k := K.carrier`・`E := K.closure` を代入し、可換性の供給元
(`𝒪_K^×` と `Ẑ`)を渡すだけである。

★核の要は `le_fixedField_commutator_iff`:

```
A ≤ (E)^{⁅Gal(E/k), Gal(E/k)⁆‾}  ↔  Gal(A/k) が可換
```

⇐ は「可換な商への写像の核は交換子群を含む」＋「`A.fixingSubgroup` は閉」
(`InfiniteGalois.fixingSubgroup_isClosed` + `Subgroup.topologicalClosure_minimal`)、
⇒ は `Γ ↠ Gal(A/k)`(`AlgEquiv.restrictNormalHom_surjective`)の核が
`A.fixingSubgroup ⊇ ⁅Γ,Γ⁆` であること。

## ★在庫から取ったもの(重複して証明していない)

* `Found/PGC/AbelianClosureSplit.lean`(B2+B3、本セッションで着地)
  * `unramifiedClosure_le_fixedField_commutator`——`K^ur ≤ K^ab`
  * `isGalois_fixedField_commutator`——`K^ab/K` は Galois
  ★B2+B3 は `abelianClosure` という名前を**意図的に置かず**
  `IntermediateField.fixedField ((commutator K.absGal).topologicalClosure)` を
  書き下していた。本ファイルはその字面をそのまま `def` にしたので、
  両者は `rfl` で繋がる(名前は 1 か所にしかない)。
* `Found/PGC/LubinTateClosure.lean::ker_restrictNormalHom_eq_fixingSubgroup`
* `Found/PGC/LubinTateClosure.lean::lubinTateClosureGalEquivUnits`(`Gal(K_π/K) ≅ 𝒪_K^×`)
* `Found/PGC/UnramifiedZhat.lean::unramifiedClosureGalEquivZHat`(`Gal(K^ur/K) ≅ Ẑ`)
  ＋ `Found/PGC/ProfiniteUnitsTorsion.lean::zhatCommGroup`(`Ẑ` は可換)
* `Found/PGC/AbelianDecomposition.lean::exists_lubinTateUnramified_decomposition`
* `Found/PGC/QpNonAbelian.lean::not_commutative_absGal`(`Γ_{ℚ_p}` は非可換)

★**`K_π` と `K^ur` の可換性は別々の在庫から来ている**——前者は `𝒪_K^×`、
後者は `Ẑ`。どちらも `mul_comm_of_mulEquiv`(可換群と同型な群は可換)を
通して同じ抽象核に入る。

★なお `K^ur ≤ K^ab` は**独立に 2 通り**得ている:
`unramifiedClosure_le_abelianClosure`(B2+B3 の交換子計算をそのまま使う)と
`unramifiedClosure_le_abelianClosure_of_zhat`(`Gal(K^ur/K) ≅ Ẑ` から抽象核へ)。
主張は同じなので後者は独立検算である。

## ★★退化の自己検査

1. **`⊥` に潰れていない**——`abelianClosure_ne_bot`。
   `K^ur ≤ K^ab` と `unramifiedClosure_ne_bot`(`Gal(K^ur/K) ≅ Ẑ` は非自明、
   `exists_zhat_ne_one`)から。
2. **`⊤` に潰れていない**——`abelianClosure_selfField_ne_top`。
   `K^ab = ⊤` なら `Γ_K` が可換になる(`forall_mul_comm_absGal_of_abelianClosure_eq_top`)が、
   `K = ℚ_p` では `not_commutative_absGal` に反する。
   ★一般の `K` について `Γ_K` が非可換であることは木にまだ無いので、
   ここは `selfField p`(= `ℚ_p`)でしか言えていない。
3. **位相的閉包を取る理由**——`⁅Γ_K,Γ_K⁆` そのもの(閉でない可能性がある)の
   固定体を取ると、無限次 Galois 対応
   (`InfiniteGalois.fixingSubgroup_fixedField`)が使えず、
   `(K^ab).fixingSubgroup = ⁅Γ_K,Γ_K⁆` が言えない。すなわち
   `Gal(K^ab/K)` の可換性(`instIsAbelianGaloisAbelianClosure`)と
   最大性(`le_abelianClosure_iff`)がどちらも壊れる。
   ★閉包を取れば `Subgroup.isClosed_topologicalClosure` が効く。
   同じ理由で `Found/PGC/TopAbelianization.lean` も `Abelianization Γ`
   (代数的商)を使っていない。
4. **`K^ab` は最大アーベル部分拡大**——`le_abelianClosure_iff` が
   両向きなので、「大きすぎ」も「小さすぎ」もしていない。

## ★逸脱の記録

なし。原典 Theorem 6.15 の記号 `K^ab` に、標準的な定義
(絶対 Galois 群の位相的交換子群の固定体)を与えただけである。
★原典は `K^ab` を「有限次アーベル拡大すべての合成体」の意味で使うが、
`le_abelianClosure_iff` はその読みと一致する
(有限次アーベル `A` は `Gal(A/K)` 可換ゆえ `A ≤ K^ab`、逆に
`K^ab` の有限次部分拡大の Galois 群は `Gal(K^ab/K)` の商で可換)。
★本ファイルは前者の形(固定体)を採った。無限次を扱うので
「合成体」を `iSup` で書くより固定体の方が Galois 対応に直結するためである。

## ★やっていないこと

* `K^ab ≤ K^LT`(B5・B6 の担当)。
* 一般の `K` についての `Γ_K` 非可換(木に無い)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-! ## 0. 抽象核——一般の Galois 拡大 `E/k` だけで閉じる部分

★本節には分岐・付値・Lubin-Tate の語彙が 1 語も出てこない。 -/

section AbstractCore

/-- 可換群と群同型な群は可換。 -/
theorem mul_comm_of_mulEquiv {G H : Type*} [Group G] [CommGroup H] (e : G ≃* H) (a b : G) :
    a * b = b * a := by
  refine e.injective ?_
  rw [map_mul, map_mul, mul_comm]

variable {k E : Type*} [Field k] [Field E] [Algebra k E]

/-- 中間体の包含に対して固定部分群は反変。 -/
theorem fixingSubgroup_le_fixingSubgroup {A B : IntermediateField k E} (h : A ≤ B) :
    B.fixingSubgroup ≤ A.fixingSubgroup := fun _σ hσ =>
  (IntermediateField.mem_fixingSubgroup_iff _ _).mpr fun x hx =>
    (IntermediateField.mem_fixingSubgroup_iff _ _).mp hσ x (h hx)

/-- `A = ⊥` なら `Gal(A/k)` は自明。 -/
theorem algEquiv_eq_one_of_eq_bot {A : IntermediateField k E} (h : A = ⊥) (σ : A ≃ₐ[k] A) :
    σ = 1 := by
  subst h
  refine AlgEquiv.ext fun x => ?_
  obtain ⟨a, ha⟩ := (IntermediateField.mem_bot).mp x.2
  have hx : x = algebraMap k (⊥ : IntermediateField k E) a := by
    apply Subtype.ext; simpa using ha.symm
  rw [hx]
  simp

/-- **抽象核(⇐ の中身)**: `Gal(A/k)` が可換なら `⁅Gal(E/k), Gal(E/k)⁆ ≤ A.fixingSubgroup`。

`A.fixingSubgroup = ker(restrictNormalHom A)` で、可換群への準同型の核は
交換子群を含むから。 -/
theorem commutator_le_fixingSubgroup_of_forall_mul_comm
    (A : IntermediateField k E) [Normal k A]
    (hab : ∀ a b : (A ≃ₐ[k] A), a * b = b * a) :
    commutator (E ≃ₐ[k] E) ≤ A.fixingSubgroup := by
  rw [← ker_restrictNormalHom_eq_fixingSubgroup A, commutator_def]
  refine Subgroup.commutator_le.mpr fun a _ b _ => ?_
  rw [MonoidHom.mem_ker, map_commutatorElement, commutatorElement_eq_one_iff_mul_comm]
  exact hab _ _

/-- **抽象核(⇒ の中身)**: `⁅Gal(E/k), Gal(E/k)⁆ ≤ A.fixingSubgroup` なら
`Gal(A/k)` は可換。

`Gal(E/k) ↠ Gal(A/k)`(`AlgEquiv.restrictNormalHom_surjective`)の核が
`A.fixingSubgroup` を含むので、交換子が 1 に落ちる。

★`lean-idioms.md` #154 に従い `restrictNormalHom` は `set` で局所定数にしている。 -/
theorem forall_mul_comm_of_commutator_le_fixingSubgroup [IsGalois k E]
    (A : IntermediateField k E) [Normal k A]
    (h : commutator (E ≃ₐ[k] E) ≤ A.fixingSubgroup) :
    ∀ a b : (A ≃ₐ[k] A), a * b = b * a := by
  set φ := AlgEquiv.restrictNormalHom (F := k) (K₁ := E) (E := (A : Type _)) with hφ
  have hker : commutator (E ≃ₐ[k] E) ≤ MonoidHom.ker φ := by
    rw [hφ, ker_restrictNormalHom_eq_fixingSubgroup]; exact h
  have hsurj : Function.Surjective φ := by
    rw [hφ]
    exact AlgEquiv.restrictNormalHom_surjective (F := k) (K₁ := (A : Type _)) E
  intro a b
  obtain ⟨α, rfl⟩ := hsurj a
  obtain ⟨β, rfl⟩ := hsurj b
  have hmem : α * β * α⁻¹ * β⁻¹ ∈ commutator (E ≃ₐ[k] E) := by
    rw [commutator_def]
    exact Subgroup.commutator_mem_commutator (Subgroup.mem_top α) (Subgroup.mem_top β)
  have h1 : φ (α * β * α⁻¹ * β⁻¹) = 1 := hker hmem
  simp only [map_mul, map_inv] at h1
  rwa [mul_inv_eq_one, mul_inv_eq_iff_eq_mul] at h1

/-- **抽象核 1**: `Gal(A/k)` が可換な正規中間体は `(E)^{⁅Γ,Γ⁆‾}` に含まれる。

`A.fixingSubgroup` は閉(`InfiniteGalois.fixingSubgroup_isClosed`)なので、
位相的閉包を取っても包含は保たれる(`Subgroup.topologicalClosure_minimal`)。 -/
theorem le_fixedField_commutator_topologicalClosure [IsGalois k E]
    (A : IntermediateField k E) [Normal k A]
    (hab : ∀ a b : (A ≃ₐ[k] A), a * b = b * a) :
    A ≤ IntermediateField.fixedField ((commutator (E ≃ₐ[k] E)).topologicalClosure) := by
  refine (IntermediateField.le_iff_le _ _).mpr ?_
  exact Subgroup.topologicalClosure_minimal _
    (commutator_le_fixingSubgroup_of_forall_mul_comm A hab)
    (InfiniteGalois.fixingSubgroup_isClosed _)

/-- **抽象核 2**: `(E)^{⁅Γ,Γ⁆‾}/k` は Galois。

`⁅Γ,Γ⁆‾` は閉かつ正規なので、無限次 Galois 対応
(`InfiniteGalois.fixingSubgroup_fixedField` と `normal_iff_isGalois`)が使える。 -/
theorem isGalois_fixedField_commutator_topologicalClosure [IsGalois k E] :
    IsGalois k (IntermediateField.fixedField ((commutator (E ≃ₐ[k] E)).topologicalClosure)) := by
  haveI : ((commutator (E ≃ₐ[k] E)).topologicalClosure).Normal :=
    Subgroup.is_normal_topologicalClosure _
  have hfix : (IntermediateField.fixedField
        ((commutator (E ≃ₐ[k] E)).topologicalClosure)).fixingSubgroup
      = (commutator (E ≃ₐ[k] E)).topologicalClosure :=
    InfiniteGalois.fixingSubgroup_fixedField
      ⟨(commutator (E ≃ₐ[k] E)).topologicalClosure, Subgroup.isClosed_topologicalClosure _⟩
  rw [← InfiniteGalois.normal_iff_isGalois, hfix]
  infer_instance

/-- **抽象核 3**: `(E)^{⁅Γ,Γ⁆‾}` に含まれる正規中間体の Galois 群は可換。

とくに `A := (E)^{⁅Γ,Γ⁆‾}` 自身に当てると「`K^ab/k` はアーベル」が出る。 -/
theorem forall_mul_comm_of_le_fixedField_commutator [IsGalois k E]
    (A : IntermediateField k E) [Normal k A]
    (hA : A ≤ IntermediateField.fixedField ((commutator (E ≃ₐ[k] E)).topologicalClosure)) :
    ∀ a b : (A ≃ₐ[k] A), a * b = b * a := by
  haveI : ((commutator (E ≃ₐ[k] E)).topologicalClosure).Normal :=
    Subgroup.is_normal_topologicalClosure _
  have hfix : (IntermediateField.fixedField
        ((commutator (E ≃ₐ[k] E)).topologicalClosure)).fixingSubgroup
      = (commutator (E ≃ₐ[k] E)).topologicalClosure :=
    InfiniteGalois.fixingSubgroup_fixedField
      ⟨(commutator (E ≃ₐ[k] E)).topologicalClosure, Subgroup.isClosed_topologicalClosure _⟩
  refine forall_mul_comm_of_commutator_le_fixingSubgroup A ?_
  refine le_trans (le_trans (Subgroup.le_topologicalClosure _) (le_of_eq hfix.symm)) ?_
  exact fixingSubgroup_le_fixingSubgroup hA

/-- ★★**抽象核 4(要)**: `(E)^{⁅Γ,Γ⁆‾}` は**最大アーベル部分拡大**である。

`k` 上正規な中間体 `A` について
`A ≤ (E)^{⁅Γ,Γ⁆‾} ↔ Gal(A/k)` が可換。 -/
theorem le_fixedField_commutator_iff [IsGalois k E]
    (A : IntermediateField k E) [Normal k A] :
    A ≤ IntermediateField.fixedField ((commutator (E ≃ₐ[k] E)).topologicalClosure) ↔
      ∀ a b : (A ≃ₐ[k] A), a * b = b * a :=
  ⟨forall_mul_comm_of_le_fixedField_commutator A,
    le_fixedField_commutator_topologicalClosure A⟩

/-- **抽象核 5(退化の自己検査)**: `(E)^{⁅Γ,Γ⁆‾} = ⊤` なら `Gal(E/k)` 自身が可換。

固定体が `⊤` なら `⁅Γ,Γ⁆‾` の元はすべての元を固定する、すなわち `⁅Γ,Γ⁆‾ = ⊥`。
★これは Galois 対応を経由しない(`mem_fixedField_iff` と `AlgEquiv.ext` だけ)。 -/
theorem forall_mul_comm_of_fixedField_commutator_eq_top
    (h : IntermediateField.fixedField ((commutator (E ≃ₐ[k] E)).topologicalClosure) = ⊤) :
    ∀ a b : (E ≃ₐ[k] E), a * b = b * a := by
  have hbot : (commutator (E ≃ₐ[k] E)).topologicalClosure = ⊥ := by
    refine le_antisymm (fun σ hσ => ?_) bot_le
    have hσ1 : σ = 1 := by
      refine AlgEquiv.ext fun x => ?_
      have hx : x ∈ IntermediateField.fixedField
          ((commutator (E ≃ₐ[k] E)).topologicalClosure) := by rw [h]; trivial
      exact (IntermediateField.mem_fixedField_iff _ _).mp hx σ hσ
    simp [hσ1]
  intro a b
  have hle : commutator (E ≃ₐ[k] E) ≤ ⊥ :=
    hbot ▸ Subgroup.le_topologicalClosure (commutator (E ≃ₐ[k] E))
  rw [commutator_def] at hle
  have hmem := Subgroup.commutator_le.mp hle a (Subgroup.mem_top a) b (Subgroup.mem_top b)
  exact commutatorElement_eq_one_iff_mul_comm.mp (by simpa using hmem)

end AbstractCore

/-! ## 1. `K^ab` の定義 -/

variable {p : ℕ} [Fact p.Prime]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`K^ab := (K̄)^{⁅Γ_K, Γ_K⁆‾}`**——`K` の最大アーベル拡大。

原典 [Yoshida08] は §6 でこの記号を説明なしに使う。木にも mathlib にも
「`K^ab` という中間体」は無かった(mathlib に在るのは
`Field.absoluteGaloisGroupAbelianization`——群の側だけ)ので、ここで立てる。

★**位相的閉包 `‾` を取る**のが要である。取らないと `(K^ab).fixingSubgroup`
が `⁅Γ_K,Γ_K⁆` に戻らず、`Gal(K^ab/K)` の可換性も最大性も言えない
(モジュール docstring「退化の自己検査 3」)。

★字面は `Found/PGC/AbelianClosureSplit.lean`(B2+B3)の書き下しと**同一**に
してあるので、両者は `rfl` で繋がる。 -/
noncomputable def abelianClosure (K : PAdicLocalField p) : IntermediateField K.carrier K.closure :=
  IntermediateField.fixedField ((commutator K.absGal).topologicalClosure)

/-- `K^ab` の定義の展開(B2+B3 の書き下しへ乗り移るため)。 -/
theorem abelianClosure_def (K : PAdicLocalField p) :
    abelianClosure K = IntermediateField.fixedField ((commutator K.absGal).topologicalClosure) :=
  rfl

/-- `K^ab/K` は(無限次)Galois。★`AbelianClosureSplit.lean` の在庫をそのまま使う。 -/
instance instIsGaloisAbelianClosure (K : PAdicLocalField p) :
    IsGalois K.carrier (abelianClosure K) :=
  isGalois_fixedField_commutator K

/-- `K^ab/K` は正規。 -/
instance instNormalAbelianClosure (K : PAdicLocalField p) :
    Normal K.carrier (abelianClosure K) :=
  IsGalois.to_normal

/-- ★**`Gal(K^ab/K)` は可換**。抽象核 3 に `A := K^ab` を当てただけ。 -/
theorem forall_mul_comm_abelianClosureGal (K : PAdicLocalField p) :
    ∀ a b : (abelianClosure K ≃ₐ[K.carrier] abelianClosure K), a * b = b * a := by
  haveI := isGalois_closure K
  exact forall_mul_comm_of_le_fixedField_commutator (abelianClosure K) le_rfl

/-- ★★`K^ab/K` は mathlib の意味で**アーベル Galois 拡大**。 -/
instance instIsAbelianGaloisAbelianClosure (K : PAdicLocalField p) :
    IsAbelianGalois K.carrier (abelianClosure K) where
  is_comm := ⟨forall_mul_comm_abelianClosureGal K⟩

/-- ★★★**`K^ab` は最大アーベル部分拡大**——`K` 上正規な `A` について
`A ≤ K^ab ↔ Gal(A/K)` が可換。

★これが「`abelianClosure` が退化していない」ことの本体である
(大きすぎれば ⇒ が、小さすぎれば ⇐ が壊れる)。 -/
theorem le_abelianClosure_iff (K : PAdicLocalField p)
    (A : IntermediateField K.carrier K.closure) [Normal K.carrier A] :
    A ≤ abelianClosure K ↔ ∀ a b : (A ≃ₐ[K.carrier] A), a * b = b * a := by
  haveI := isGalois_closure K
  exact le_fixedField_commutator_iff A

/-- 可換な `Gal(A/K)` を持つ正規中間体は `K^ab` に含まれる(⇐ 向きの単体形)。 -/
theorem le_abelianClosure (K : PAdicLocalField p)
    (A : IntermediateField K.carrier K.closure) [Normal K.carrier A]
    (hab : ∀ a b : (A ≃ₐ[K.carrier] A), a * b = b * a) :
    A ≤ abelianClosure K :=
  (le_abelianClosure_iff K A).mpr hab

/-! ## 2. `K^ur ≤ K^ab` -/

/-- ★**`K^ur ≤ K^ab`**。`Found/PGC/AbelianClosureSplit.lean`(B2+B3)の
`unramifiedClosure_le_fixedField_commutator` に `rfl` で乗り移るだけ。 -/
theorem unramifiedClosure_le_abelianClosure (K : PAdicLocalField p) :
    unramifiedClosure K ≤ abelianClosure K :=
  unramifiedClosure_le_fixedField_commutator K

/-- ★**`Gal(K^ur/K)` は可換**——`Gal(K^ur/K) ≅ Ẑ`(`unramifiedClosureGalEquivZHat`)と
`Ẑ` の可換性(`zhatCommGroup`)から。

★これは `K_π` 側(`𝒪_K^×`)とは**別の在庫**である。 -/
theorem forall_mul_comm_unramifiedClosureGal (K : PAdicLocalField p) :
    ∀ a b : (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)), a * b = b * a :=
  fun a b => mul_comm_of_mulEquiv (unramifiedClosureGalEquivZHat K).toMulEquiv a b

/-- ★`K^ur ≤ K^ab` の**独立な第 2 証明**——`Ẑ` の可換性から抽象核へ。

`unramifiedClosure_le_abelianClosure`(B2+B3 の交換子計算)と主張は同じで、
証明の経路だけが違う。★片方が壊れてももう片方が残る。 -/
theorem unramifiedClosure_le_abelianClosure_of_zhat (K : PAdicLocalField p) :
    unramifiedClosure K ≤ abelianClosure K :=
  le_abelianClosure K (unramifiedClosure K) (forall_mul_comm_unramifiedClosureGal K)

/-! ## 3. `K_π ≤ K^ab` -/

section LubinTate

variable (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))

/-- ★**`Gal(K_π/K)` は可換**——`Gal(K_π/K) ≅ 𝒪_K^×`
(`lubinTateClosureGalEquivUnits`)と `Units` が `CommGroup` であることから。

★これは `K^ur` 側(`Ẑ`)とは**別の在庫**である。 -/
theorem forall_mul_comm_lubinTateClosureGal :
    ∀ a b : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ≃ₐ[K.carrier]
      lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf), a * b = b * a :=
  fun a b => mul_comm_of_mulEquiv (lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf) a b

/-- ★★**`K_π ≤ K^ab`**。 -/
theorem lubinTateClosure_le_abelianClosure :
    lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ≤ abelianClosure K := by
  haveI := normal_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
  exact le_abelianClosure K _ (forall_mul_comm_lubinTateClosureGal K hq hπmax hπne0 f hf0 hf1 hf)

/-- ★★★**`K^LT = K_π · K^ur ≤ K^ab`**。

★`≤` であって `=` ではない(逆向きは B5・B6)。 -/
theorem lubinTateClosure_sup_unramifiedClosure_le_abelianClosure :
    (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
        IntermediateField K.carrier K.closure) ≤ abelianClosure K :=
  sup_le (lubinTateClosure_le_abelianClosure K hq hπmax hπne0 f hf0 hf1 hf)
    (unramifiedClosure_le_abelianClosure K)

end LubinTate

/-! ## 4. 仮説なしの形 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`K^LT ≤ K^ab`(仮説なし)**——[Yoshida08] Theorem 6.15 の**半分**。

任意の p 進局所体 `K` について、`K̄` の中に `K` 上正規な `K_π` が在って

* `K_π ⊓ K^ur = K`
* `Gal(K_π/K) ≅ 𝒪_K^×`
* **`K_π · K^ur ≤ K^ab`**

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

★★**本ノード(B1)は `≤` しか主張しない。** 原典の等号 `K^LT = K^ab` は
逆向き `K^ab ≤ K^LT`(B5: `E_σ ⊆ K^ram_x`、B6: 合成)が入って初めて閉じる。
★取り違えないこと。

★`Λ_∞` の選び方(`π`・`f`)に依らない形にするため、`K_π` は `∃` の内側に
閉じ込めてある(`exists_lubinTateUnramified_decomposition` と同じ流儀)。 -/
theorem exists_lubinTate_le_abelianClosure (K : PAdicLocalField p) :
    ∃ E : IntermediateField K.carrier K.closure,
      Normal K.carrier E ∧
      E ⊓ unramifiedClosure K = ⊥ ∧
      Nonempty ((E ≃ₐ[K.carrier] E) ≃* (𝒪[K.carrier])ˣ) ∧
      (E ⊔ unramifiedClosure K : IntermediateField K.carrier K.closure) ≤ abelianClosure K := by
  obtain ⟨E, hEn, hEinf, ⟨e⟩, hprod⟩ := exists_lubinTateUnramified_decomposition K
  refine ⟨E, hEn, hEinf, ⟨e⟩, sup_le ?_ (unramifiedClosure_le_abelianClosure K)⟩
  haveI := hEn
  exact le_abelianClosure K E (fun a b => mul_comm_of_mulEquiv e a b)

def exists_lubinTate_le_abelianClosure.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-! ## 5. 退化の自己検査 -/

/-- `K^ur ≠ K`——`Gal(K^ur/K) ≅ Ẑ` は非自明(`exists_zhat_ne_one`)だから。 -/
theorem unramifiedClosure_ne_bot (K : PAdicLocalField p) : unramifiedClosure K ≠ ⊥ := by
  intro h
  obtain ⟨x, hx⟩ := exists_zhat_ne_one
  exact hx (by
    have := algEquiv_eq_one_of_eq_bot h ((unramifiedClosureGalEquivZHat K).toMulEquiv.symm x)
    simpa using congrArg (unramifiedClosureGalEquivZHat K).toMulEquiv this)

/-- ★**`K^ab ≠ K`**(`⊥` に潰れていない)——`K ⊊ K^ur ≤ K^ab`。 -/
theorem abelianClosure_ne_bot (K : PAdicLocalField p) : abelianClosure K ≠ ⊥ := by
  intro h
  exact unramifiedClosure_ne_bot K
    (le_antisymm (h ▸ unramifiedClosure_le_abelianClosure K) bot_le)

/-- `K^ab = K̄` なら `Γ_K` は可換(抽象核 5)。 -/
theorem forall_mul_comm_absGal_of_abelianClosure_eq_top (K : PAdicLocalField p)
    (h : abelianClosure K = ⊤) : ∀ a b : K.absGal, a * b = b * a :=
  forall_mul_comm_of_fixedField_commutator_eq_top h

/-- ★**`K^ab ≠ K̄`**(`⊤` にも潰れていない)——`K = ℚ_p` の場合。

`Γ_{ℚ_p}` は非可換(`Found/PGC/QpNonAbelian.lean::not_commutative_absGal`)。
★一般の `K` について `Γ_K` 非可換は木にまだ無いので、ここは `ℚ_p` に限る。 -/
theorem abelianClosure_selfField_ne_top (p : ℕ) [Fact p.Prime] :
    abelianClosure (selfField p) ≠ ⊤ := fun h =>
  not_commutative_absGal p (forall_mul_comm_absGal_of_abelianClosure_eq_top (selfField p) h)

end ABC3.Found.PGC
