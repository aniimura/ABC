import ABC3.Found.PGC.Theorem42Bijectivity
import ABC3.Found.PGC.SubgroupCorrespondenceConstruction

/-!
# [pGC] Theorem 4.2 単射性の要 —— 中心化群の条件を「有限次の 1 文」へ落とす

`Found/PGC/Theorem42Bijectivity.lean` は pGC Theorem 4.2 の単射性が

  `CentralizerActsTriviallyOnBase K`  （= `C_{Γ_Qp}(Γ_K) ⊆ Γ_K` の弱形）

と **同値** であることを示した（`injective_naturalOuterIso_iff`）。本ファイルは
その 1 文を**代数閉包を 1 度も見ない有限次の条件**に落とす。

原文 (pGC 物理 p.7, Theorem 4.2 の証明の括弧書き):

> (Alternatively, the more group-theoretically
> oriented reader may prefer to regard the injectivity of this morphism as a consequence of
> the fact that the centralizer of ΓK in ΓQp is trivial.)

## 何を得たか（測ってあること・測っていないこと）

| 段 | 結果 |
|---|---|
| 抽象核 | `fixed_of_comm_of_fixed` —— 「族と可換な写像は共通不動点を保つ」。★型クラスが 1 つも要らない。 |
| Galois 層 | `apply_mem_intermediateField_of_comm` —— Galois 拡大 `E/F` で `Γ_F` と可換な `σ` は**すべての中間体を保つ**。 |
| 障害 | `exists_root_adjoin_of_centralizing` —— 中心化する `σ` があれば、`q(y)=0` なる任意の `y` について `ρ(q)` が `K(y)` の中に根を持つ。 |
| 判定条件 | `centralizerActsTriviallyOnBase_of_minpoly_witness` / `_of_root_witness` / `_of_sq_witness`。 |
| 有限次還元 | `centralizerActsTriviallyOnBase_of_finiteLevel` —— 有限次正規部分拡大 `L/K` 1 つで足りる。 |
| 同上（群の言葉） | `centralizerActsTriviallyOnBase_of_finiteLevelAlg` —— `L/ℚ_p` を Galois に取れば条件は `C_{Gal(L/ℚ_p)}(Gal(L/K)) ⊆ Gal(L/K)`。★有限群 1 個の計算。 |

★**未証明のまま残ったもの**: `CentralizerActsTriviallyOnBase K` そのもの。
本ファイルは**十分条件を 4 通り**与えただけであり、どの十分条件も特定の `K` について
「そういう `y`（または `L`）が実在する」ことを言うには p 進の入力が要る。
★`sorry` は 1 つも無い。

## ★なぜ「形式的な証明」は存在し得ないか（`Check/PGC/CentralizerNeedsInput.lean`）

同じ文を「基礎体 `k`・その有限次拡大 `K`・`K` の代数閉包」の言葉に一般化すると**偽**である
（`k = ℝ`, `K = ℂ` が反例）。すなわち本条件は p 進性を本質的に使う。
★この観測は `Check/PGC/CentralizerNeedsInput.lean` に反例つきで置いた。

## ★原典より短い道（見つけたもの）

原典（および素朴な期待）は「中心化群が自明」を **`Z(Γ_F) = 1` 型の入力**から出す。
本ファイルが使うのは **無限 Galois 対応 1 本だけ**である:

  σ が Γ_K と可換 ⟹ σ は `Γ_L`(`L ⊇ K`) を保つ ⟹ σ(L) ⊆ L（`L = fixedField Γ_L` だから）。

これで `σ y ∈ K(y)` が**任意の** `y` について出る。この 1 行が本ファイルの全部を支えている。
★分岐も付値も Frobenius も 1 度も出てこない。

## ★測っていないこと（正直に）

* 具体的な `K`（`K ≠ ℚ_p`）について判定条件の witness を**構成していない**。
  ゆえに `Theorem42Bijectivity.lean` の非空虚性 witness（`K = ℚ_p`）は**弱いまま**である。
* `_of_sq_witness` は**万能ではない**: 例えば `K/ℚ_p` が不分岐 2 次のとき、Frobenius は
  `K^×/(K^×)^2` に自明に働く（剰余体で 2 乗と Frobenius が可換だから）ので、
  平方根だけでは `ρ ≠ 1` を排除できない。一般次数版 `_of_root_witness` が要る理由である。
  ★この観察自体は本リポジトリでは**形式化していない**（紙の上の注意）。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC IntermediateField

/-! ## §1 抽象核 —— 型クラスが 1 つも要らない 1 文

★ここには群も環も体も出てこない。「写像の族と可換な写像は、その族の共通不動点を保つ」。 -/

section AbstractCore

variable {X : Type*} {ι : Type*}

/-- ★★★**抽象核**。写像の族 `G : ι → (X → X)` と可換な `σ : X → X` は、
族の任意の部分族 `s` の共通不動点を保つ。

★型クラスの仮定が **0 個**である。分岐・付値・Galois どころか、群も環も出てこない。 -/
theorem fixed_of_comm_of_fixed (G : ι → (X → X)) (σ : X → X)
    (hσ : ∀ (i : ι) (x : X), σ (G i x) = G i (σ x)) {y : X} {s : Set ι}
    (hy : ∀ i ∈ s, G i y = y) : ∀ i ∈ s, G i (σ y) = σ y := by
  intro i hi
  rw [← hσ i y, hy i hi]

/-- 抽象核のモノイド作用版。 -/
theorem smul_fixed_of_comm_of_fixed {M : Type*} [Monoid M] [MulAction M X] (σ : X → X)
    (hσ : ∀ (g : M) (x : X), σ (g • x) = g • σ x) {y : X} {s : Set M}
    (hy : ∀ g ∈ s, g • y = y) : ∀ g ∈ s, g • σ y = σ y :=
  fixed_of_comm_of_fixed (fun g : M => fun x : X => g • x) σ hσ hy

end AbstractCore

/-! ## §2 Galois 層 —— p 進も分岐も出てこない

`E/F` を（無限次でよい）Galois 拡大とする。`Γ_F = Aut(E/F)` と可換な写像は
**すべての中間体を保つ**。証明は無限 Galois 対応（`fixedField (fixingSubgroup L) = L`）
1 本だけである。 -/

section GaloisLayer

variable {F E : Type*} [Field F] [Field E] [Algebra F E] [IsGalois F E]

/-- ★★`Aut(E/F)` の全元と可換な `σ : E → E` は、任意の中間体 `L` を保つ。

★`σ` に環準同型であることすら要求していない（単なる写像でよい）。 -/
theorem apply_mem_intermediateField_of_comm (σ : E → E)
    (hσ : ∀ (g : E ≃ₐ[F] E) (x : E), σ (g x) = g (σ x))
    (L : IntermediateField F E) {y : E} (hy : y ∈ L) : σ y ∈ L := by
  have h : ∀ g ∈ L.fixingSubgroup, g (σ y) = σ y := by
    refine fixed_of_comm_of_fixed (fun g : E ≃ₐ[F] E => fun x : E => g x) σ hσ ?_
    intro g hg
    exact (IntermediateField.mem_fixingSubgroup_iff L g).mp hg y hy
  rw [← InfiniteGalois.fixedField_fixingSubgroup L]
  simpa using h

/-- 単生成の中間体への特殊化: `σ y ∈ F(y)`。 -/
theorem apply_mem_adjoin_simple_of_comm (σ : E → E)
    (hσ : ∀ (g : E ≃ₐ[F] E) (x : E), σ (g x) = g (σ x)) (y : E) :
    σ y ∈ IntermediateField.adjoin F ({y} : Set E) :=
  apply_mem_intermediateField_of_comm σ hσ _ (mem_adjoin_simple_self F y)

end GaloisLayer

variable {p : ℕ} [Fact p.Prime]

/-! ## §3 pGC 層 -/

/-- `Γ_K` を中心化する `σ` は `K̄/K` のすべての中間体を保つ。 -/
theorem centralizing_apply_mem_intermediateField {K : PAdicLocalField p}
    (σ : K.closure → K.closure)
    (hcomm : ∀ (g : K.absGal) (x : K.closure), σ (g x) = g (σ x))
    (L : IntermediateField K.carrier K.closure) {y : K.closure} (hy : y ∈ L) : σ y ∈ L := by
  haveI := isGalois_closure K
  exact apply_mem_intermediateField_of_comm σ hcomm L hy

/-- `σ` が `Γ_K` と可換なら `σ⁻¹` もそうである。 -/
theorem centralizing_symm_comm {K : PAdicLocalField p} (σ : K.closure ≃+* K.closure)
    (hcomm : ∀ (g : K.absGal) (x : K.closure), σ (g x) = g (σ x)) :
    ∀ (g : K.absGal) (x : K.closure), σ.symm (g x) = g (σ.symm x) := by
  intro g x
  have h := hcomm g (σ.symm x)
  rw [RingEquiv.apply_symm_apply] at h
  rw [← h, RingEquiv.symm_apply_apply]

/-! ## §4 障害 —— 中心化する `σ` が在れば何が強制されるか -/

/-- `σ` が `K` の上で `ρ` として働くなら、`q(y) = 0` から `ρ(q)(σ y) = 0` が出る。 -/
theorem centralizing_aeval_map_eq_zero {K : PAdicLocalField p}
    {ρ : K.carrier ≃ₐ[ℚ_[p]] K.carrier} {σ : K.closure ≃+* K.closure}
    (hfwd : ∀ x : K.carrier, σ (algebraMap K.carrier K.closure x)
      = algebraMap K.carrier K.closure (ρ x))
    {y : K.closure} {q : Polynomial K.carrier} (hq : (Polynomial.aeval y q : K.closure) = 0) :
    (Polynomial.aeval (σ y) (q.map ρ.toRingEquiv.toRingHom) : K.closure) = 0 := by
  have hcompose : σ.toRingHom.comp (algebraMap K.carrier K.closure)
      = (algebraMap K.carrier K.closure).comp ρ.toRingEquiv.toRingHom :=
    RingHom.ext hfwd
  have h1 := Polynomial.hom_eval₂ q (algebraMap K.carrier K.closure) σ.toRingHom y
  rw [← Polynomial.aeval_def, hq, map_zero, hcompose] at h1
  rw [Polynomial.aeval_def, Polynomial.eval₂_map]
  exact h1.symm

/-- ★★★**障害（本ファイルの数学的な中身）**。

`Γ_K` を中心化し `K` の上で `ρ` として働く `σ` が在るとする。このとき
`q(y) = 0` なる任意の `y ∈ K̄` と `q ∈ K[X]` について、
`ρ` で係数をひねった `ρ(q)` は **`K(y)` の中に根を持つ**。

★根は `σ y` そのものである。★`K(y)` の外に出られないことが要点で、
これは「`σ` が中間体を保つ」（§2）から出る。 -/
theorem exists_root_adjoin_of_centralizing {K : PAdicLocalField p}
    {ρ : K.carrier ≃ₐ[ℚ_[p]] K.carrier} {σ : K.closure ≃+* K.closure}
    (hfwd : ∀ x : K.carrier, σ (algebraMap K.carrier K.closure x)
      = algebraMap K.carrier K.closure (ρ x))
    (hcomm : ∀ (g : K.absGal) (x : K.closure), σ (g x) = g (σ x))
    (y : K.closure) {q : Polynomial K.carrier} (hq : (Polynomial.aeval y q : K.closure) = 0) :
    ∃ z ∈ IntermediateField.adjoin K.carrier ({y} : Set K.closure),
      (Polynomial.aeval z (q.map ρ.toRingEquiv.toRingHom) : K.closure) = 0 :=
  ⟨σ y,
    centralizing_apply_mem_intermediateField σ hcomm _
      (mem_adjoin_simple_self K.carrier y),
    centralizing_aeval_map_eq_zero hfwd hq⟩

def exists_root_adjoin_of_centralizing.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-! ## §5 判定条件 —— 代数閉包を見ずに `CentralizerActsTriviallyOnBase K` を出す -/

/-- ★★★**多項式による判定条件**。`ρ ≠ id` のたびに
「`q(y) = 0` なのに `ρ(q)` は `K(y)` に根を持たない」ような `(y, q)` が 1 組見つかれば、
`CentralizerActsTriviallyOnBase K` が成り立つ。

★これは §4 の障害の対偶である。 -/
theorem centralizerActsTriviallyOnBase_of_root_witness (K : PAdicLocalField p)
    (h : ∀ ρ : K.carrier ≃ₐ[ℚ_[p]] K.carrier, ρ ≠ AlgEquiv.refl →
      ∃ (y : K.closure) (q : Polynomial K.carrier),
        (Polynomial.aeval y q : K.closure) = 0 ∧
        ∀ z ∈ IntermediateField.adjoin K.carrier ({y} : Set K.closure),
          (Polynomial.aeval z (q.map ρ.toRingEquiv.toRingHom) : K.closure) ≠ 0) :
    CentralizerActsTriviallyOnBase K := by
  intro ρ σ hfwd hcomm
  by_contra hne
  obtain ⟨y, q, hq0, hno⟩ := h ρ hne
  obtain ⟨z, hz, hz0⟩ := exists_root_adjoin_of_centralizing hfwd hcomm y hq0
  exact hno z hz hz0

def centralizerActsTriviallyOnBase_of_root_witness.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-- ★**最小多項式版**（実際に使う形）。`ρ ≠ id` のたびに、`y ∈ K̄` で
「`y` の `K` 上の最小多項式を `ρ` でひねったものが `K(y)` に根を持たない」ものが在ればよい。 -/
theorem centralizerActsTriviallyOnBase_of_minpoly_witness (K : PAdicLocalField p)
    (h : ∀ ρ : K.carrier ≃ₐ[ℚ_[p]] K.carrier, ρ ≠ AlgEquiv.refl →
      ∃ y : K.closure, ∀ z ∈ IntermediateField.adjoin K.carrier ({y} : Set K.closure),
        (Polynomial.aeval z ((minpoly K.carrier y).map ρ.toRingEquiv.toRingHom) : K.closure)
          ≠ 0) :
    CentralizerActsTriviallyOnBase K := by
  refine centralizerActsTriviallyOnBase_of_root_witness K fun ρ hρ => ?_
  obtain ⟨y, hy⟩ := h ρ hρ
  exact ⟨y, minpoly K.carrier y, minpoly.aeval K.carrier y, hy⟩

/-- ★**平方根版**（いちばん手で確かめやすい形）。`ρ ≠ id` のたびに `a ∈ K` と
`y ∈ K̄`（`y² = a`）で「`ρ(a)` が `K(√a)` の中で平方でない」ものが在ればよい。

☆★**万能ではない**: `K/ℚ_p` が不分岐 2 次のとき Frobenius は `K^×/(K^×)²` に自明に働くので、
この形では排除できない。一般次数の `_of_root_witness` が要る。 -/
theorem centralizerActsTriviallyOnBase_of_sq_witness (K : PAdicLocalField p)
    (h : ∀ ρ : K.carrier ≃ₐ[ℚ_[p]] K.carrier, ρ ≠ AlgEquiv.refl →
      ∃ (a : K.carrier) (y : K.closure), y ^ 2 = algebraMap K.carrier K.closure a ∧
        ∀ z ∈ IntermediateField.adjoin K.carrier ({y} : Set K.closure),
          z ^ 2 ≠ algebraMap K.carrier K.closure (ρ a)) :
    CentralizerActsTriviallyOnBase K := by
  intro ρ σ hfwd hcomm
  by_contra hne
  obtain ⟨a, y, hy, hno⟩ := h ρ hne
  have hmem := centralizing_apply_mem_intermediateField (K := K) σ hcomm
    (IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (mem_adjoin_simple_self K.carrier y)
  refine hno _ hmem ?_
  rw [← map_pow σ y 2, hy, hfwd a]

/-! ## §6 有限次への還元 —— 代数閉包が消える -/

/-- `Γ_K` を中心化する `σ` の、中間体 `L` への制限。§2 で `σ` も `σ⁻¹` も `L` を保つので
`L ≃+* L` が立つ。 -/
noncomputable def centralizingRestrict {K : PAdicLocalField p} (σ : K.closure ≃+* K.closure)
    (hcomm : ∀ (g : K.absGal) (x : K.closure), σ (g x) = g (σ x))
    (L : IntermediateField K.carrier K.closure) : L ≃+* L where
  toFun y := ⟨σ y, centralizing_apply_mem_intermediateField σ hcomm L y.2⟩
  invFun y := ⟨σ.symm y, centralizing_apply_mem_intermediateField σ.symm
    (centralizing_symm_comm σ hcomm) L y.2⟩
  left_inv y := by ext; simp
  right_inv y := by ext; simp
  map_mul' a b := by ext; simp
  map_add' a b := by ext; simp

@[simp] theorem centralizingRestrict_coe {K : PAdicLocalField p} (σ : K.closure ≃+* K.closure)
    (hcomm : ∀ (g : K.absGal) (x : K.closure), σ (g x) = g (σ x))
    (L : IntermediateField K.carrier K.closure) (y : L) :
    ((centralizingRestrict σ hcomm L y : L) : K.closure) = σ y := rfl

/-- 制限は `Gal(L/K)` を中心化する（`L/K` が正規なら）。

★`Γ_K ↠ Gal(L/K)`（`AlgEquiv.restrictNormalHom_surjective`）で持ち上げて `hcomm` を使う。 -/
theorem centralizingRestrict_comm {K : PAdicLocalField p} (σ : K.closure ≃+* K.closure)
    (hcomm : ∀ (g : K.absGal) (x : K.closure), σ (g x) = g (σ x))
    (L : IntermediateField K.carrier K.closure) [Normal K.carrier L]
    (g' : L ≃ₐ[K.carrier] L) (y : L) :
    centralizingRestrict σ hcomm L (g' y) = g' (centralizingRestrict σ hcomm L y) := by
  obtain ⟨g, hg⟩ := AlgEquiv.restrictNormalHom_surjective (F := K.carrier) (K₁ := L)
    (E := K.closure) g'
  apply Subtype.ext
  have hres : ∀ z : L, ((g' z : L) : K.closure) = g (z : K.closure) := by
    intro z
    rw [← hg]
    exact AlgEquiv.restrictNormal_commutes g L z
  show (σ ((g' y : L) : K.closure)) = ((g' (centralizingRestrict σ hcomm L y) : L) : K.closure)
  rw [hres y, hres (centralizingRestrict σ hcomm L y), centralizingRestrict_coe]
  exact hcomm g y

/-- 制限は `K` の上で `ρ` として働く。 -/
theorem centralizingRestrict_algebraMap {K : PAdicLocalField p}
    {ρ : K.carrier ≃ₐ[ℚ_[p]] K.carrier} (σ : K.closure ≃+* K.closure)
    (hcomm : ∀ (g : K.absGal) (x : K.closure), σ (g x) = g (σ x))
    (hfwd : ∀ x : K.carrier, σ (algebraMap K.carrier K.closure x)
      = algebraMap K.carrier K.closure (ρ x))
    (L : IntermediateField K.carrier K.closure) (x : K.carrier) :
    centralizingRestrict σ hcomm L (algebraMap K.carrier L x)
      = algebraMap K.carrier L (ρ x) := by
  apply Subtype.ext
  rw [centralizingRestrict_coe]
  simpa using hfwd x

/-- ★★★★**有限次への還元**。

`ρ ≠ id` のたびに、`K̄/K` の**正規な**中間体 `L` を 1 つ選んで

  「`K` の上で `ρ` として働き、かつ `Gal(L/K)` を中心化する `s : L ≃+* L` は存在しない」

が言えれば、`CentralizerActsTriviallyOnBase K` が成り立つ。

★★**使い方**: `L` を `ℚ_p` 上有限次 Galois に取れば、`s` は `Gal(L/ℚ_p) =: G` の元であり、
条件は `H := Gal(L/K)` について `C_G(H) ⊆ H`（`s ∉ H` なのに `s ∈ C_G(H)` は無い）——
すなわち **有限群 1 個の中心化群の計算**になる。代数閉包も無限 Galois 群も出てこない。

★`L` の有限次性は証明には要らないので仮定していない（`L = K̄` を入れれば元の文に戻るだけ）。 -/
theorem centralizerActsTriviallyOnBase_of_finiteLevel (K : PAdicLocalField p)
    (h : ∀ ρ : K.carrier ≃ₐ[ℚ_[p]] K.carrier, ρ ≠ AlgEquiv.refl →
      ∃ (L : IntermediateField K.carrier K.closure) (_ : Normal K.carrier L),
        ∀ s : L ≃+* L,
          (∀ x : K.carrier, s (algebraMap K.carrier L x) = algebraMap K.carrier L (ρ x)) →
          (∀ (g : L ≃ₐ[K.carrier] L) (y : L), s (g y) = g (s y)) → False) :
    CentralizerActsTriviallyOnBase K := by
  intro ρ σ hfwd hcomm
  by_contra hne
  obtain ⟨L, hnorm, hno⟩ := h ρ hne
  haveI := hnorm
  exact hno (centralizingRestrict σ hcomm L)
    (centralizingRestrict_algebraMap σ hcomm hfwd L)
    (centralizingRestrict_comm σ hcomm L)

def centralizerActsTriviallyOnBase_of_finiteLevel.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-- ★★★**有限次還元の `Gal(L/ℚ_p)` 版**（数論の言葉に直したもの）。

`s` を `L ≃ₐ[ℚ_p] L`、すなわち `L/ℚ_p` が Galois なら `G := Gal(L/ℚ_p)` の元として
取ってよい。条件は `H := Gal(L/K)` について

  「`C_G(H)` の元で `K` の上に `ρ (≠ id)` を誘導するものは無い」

——`H` の元は `K` の上で恒等だから、これは **`C_G(H) ⊆ H`** の言い換えである。
★すなわち `CentralizerActsTriviallyOnBase K` は
**有限群 1 個の中心化群の計算**に落ちる。

★`s` が `ℚ_p` を固定することは仮定ではなく**結論**である（`ρ` が `ℚ_p`-代数同型なので
`s` は `ℚ_p` 上で自動的に恒等）。ゆえにこの版は §6 の `L ≃+* L` 版と**同じ強さ**を持つ。 -/
theorem centralizerActsTriviallyOnBase_of_finiteLevelAlg (K : PAdicLocalField p)
    (h : ∀ ρ : K.carrier ≃ₐ[ℚ_[p]] K.carrier, ρ ≠ AlgEquiv.refl →
      ∃ (L : IntermediateField K.carrier K.closure) (_ : Normal K.carrier L),
        ∀ s : L ≃ₐ[ℚ_[p]] L,
          (∀ x : K.carrier, s (algebraMap K.carrier L x) = algebraMap K.carrier L (ρ x)) →
          (∀ (g : L ≃ₐ[K.carrier] L) (y : L), s (g y) = g (s y)) → False) :
    CentralizerActsTriviallyOnBase K := by
  refine centralizerActsTriviallyOnBase_of_finiteLevel K fun ρ hρ => ?_
  obtain ⟨L, hnorm, hno⟩ := h ρ hρ
  refine ⟨L, hnorm, fun s hs hcomm => ?_⟩
  have hcommutes : ∀ c : ℚ_[p], s (algebraMap ℚ_[p] L c) = algebraMap ℚ_[p] L c := by
    intro c
    rw [IsScalarTower.algebraMap_apply ℚ_[p] K.carrier L c, hs, ρ.commutes c]
  exact hno (AlgEquiv.ofRingEquiv hcommutes) hs hcomm

def centralizerActsTriviallyOnBase_of_finiteLevelAlg.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-! ## §7 古典的な道（`Z(Γ_F) = 1` 経由）の**純群論部分**

原典の括弧書きを字面通りに追う古典的な道は次の 4 歩である。

1. `σ` が `Γ_K` を中心化する ⟹ `σ` は `Γ_K` を正規化する ⟹ `σ(K) = K`、
   すなわち `ρ ∈ Aut(K/ℚ_p)`（★本ファイル §2 が `σ(L) = L` を全中間体について与えるので、
   この歩は済んでいる）。
2. `F := K^ρ`（`ρ` の不動体）とおくと `K/F` は巡回 Galois で `Gal(K/F) = ⟨ρ⟩`。
   `Γ_F` は `Γ_K` と `σ` で**群として**生成される（`τ|_K = ρ^k` から `τ σ^{-k} ∈ Γ_K`）。
3. ゆえに `σ ∈ Z(Γ_F)`（★この歩が下の `mem_center_of_centralizes_of_zpow_generates`。
   純群論であり、体も付値も出てこない）。
4. `F` が p 進局所体なら `Z(Γ_F) = 1` ⟹ `σ = 1` ⟹ `ρ = id`。

★**4 が研究レベルの入力**である（p 進局所体の絶対 Galois 群の中心が自明であること）。
★**2 も本リポジトリには無い**（`K^ρ` を `PAdicLocalField` として立て、`Γ_F = Γ_K⟨σ⟩` を
無限 Galois 対応で示す必要がある）。本ファイルが与えるのは 1 と 3 だけである。
★§6 の有限次還元は **2・4 を回避する別の道**であり、そちらは「有限群 1 個の中心化群」で済む。 -/

/-- ★★**古典的な道の 3 歩目（純群論）**。`σ` が `H` を中心化し、
`G` が `H` と `σ` の冪で尽くされるなら `σ ∈ Z(G)`。

★体・付値・Galois の語彙が 1 語も出てこない。証明は `Commute` の 2 行。 -/
theorem mem_center_of_centralizes_of_zpow_generates {G : Type*} [Group G] {H : Subgroup G}
    {σ : G} (hcent : ∀ h ∈ H, σ * h = h * σ)
    (hgen : ∀ g : G, ∃ h ∈ H, ∃ n : ℤ, g = h * σ ^ n) :
    σ ∈ Subgroup.center G := by
  rw [Subgroup.mem_center_iff]
  intro g
  obtain ⟨h, hh, n, rfl⟩ := hgen g
  have h1 : Commute σ h := hcent h hh
  have h2 : Commute σ (σ ^ n) := (Commute.refl σ).zpow_right n
  exact ((h1.mul_right h2).eq).symm

/-! ## §8 単射性への直結 -/

/-- ★有限次の条件から、実物の分岐濾過での単射性が直接出る
（`Theorem42Bijectivity.lean::injective_naturalOuterIso_ramificationFiltration` と合成）。 -/
theorem injective_naturalOuterIso_of_finiteLevel {K K' : PAdicLocalField p}
    (h : ∀ ρ : K.carrier ≃ₐ[ℚ_[p]] K.carrier, ρ ≠ AlgEquiv.refl →
      ∃ (L : IntermediateField K.carrier K.closure) (_ : Normal K.carrier L),
        ∀ s : L ≃+* L,
          (∀ x : K.carrier, s (algebraMap K.carrier L x) = algebraMap K.carrier L (ρ x)) →
          (∀ (g : L ≃ₐ[K.carrier] L) (y : L), s (g y) = g (s y)) → False) :
    Function.Injective (naturalOuterIso (ramificationFiltration p)
      (isNaturalFiltration_ramificationFiltration p) (K := K) (K' := K')) :=
  injective_naturalOuterIso_ramificationFiltration
    (centralizerActsTriviallyOnBase_of_finiteLevel K h)

def injective_naturalOuterIso_of_finiteLevel.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

#print axioms fixed_of_comm_of_fixed
#print axioms smul_fixed_of_comm_of_fixed
#print axioms apply_mem_intermediateField_of_comm
#print axioms apply_mem_adjoin_simple_of_comm
#print axioms centralizing_apply_mem_intermediateField
#print axioms centralizing_symm_comm
#print axioms centralizing_aeval_map_eq_zero
#print axioms exists_root_adjoin_of_centralizing
#print axioms centralizerActsTriviallyOnBase_of_root_witness
#print axioms centralizerActsTriviallyOnBase_of_minpoly_witness
#print axioms centralizerActsTriviallyOnBase_of_sq_witness
#print axioms centralizingRestrict
#print axioms centralizingRestrict_coe
#print axioms centralizingRestrict_comm
#print axioms centralizingRestrict_algebraMap
#print axioms centralizerActsTriviallyOnBase_of_finiteLevel
#print axioms centralizerActsTriviallyOnBase_of_finiteLevelAlg
#print axioms mem_center_of_centralizes_of_zpow_generates
#print axioms injective_naturalOuterIso_of_finiteLevel

end ABC3.Found.PGC
