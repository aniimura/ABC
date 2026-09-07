import ABC3.Found.PGC.LubinTateCompletionDegree
import ABC3.Found.PGC.LubinTateReciprocityIsomorphism

/-!
# `Gal(K̂^m_f/K̂^{ur}) ≃ (𝒪_K/π^n)^×` —— Yoshida 2008 Proposition 4.4(iii) の `L = K̂` の場合

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Proposition 4.4(物理 p.7)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-4.html` の `#prop-4-4`。

原文 (Yoshida08 p.7):
> Proposition 4.4. Let m ≥ 1 and f ∈ O[scr]_L[X] as above, with the linear coefficient π. (i) The set µ_f,m is an O[scr]-module by +_F_f and [·]_f. For any α ∈ µ^×_f,m := µ_f,m \ µ_f,m−1, the following is an isomorphism of O[scr]-modules: O[scr]/p[frak]^m ∋ a mod p[frak]^m −→ [a]_f(α) ∈ µ_f,m. (ii) If α ∈ µ^×_f,m, then L^m_f = L(α), N_L^m_f/L(−α) = π^ϕ^m−1 and α is a uniformizer of L^m_f. The L^m_f/L is totally ramified Galois extension of degree |µ^×_f,m| = q^m−1(q − 1). (iii) We have canonical isomorphisms of abelian groups: ρ_f,m : Gal(L^m_f/L) ∼ =−→ Aut_O[scr](µ_f,m) ∼ =−→ (O[scr]/p[frak]^m)^×. (α → [u]_f(α), ∀α ∈ µ_f,m) −→ u mod p[frak]^m

★★本ノードは Proposition 4.4(iii) の **`L = K̂`(＝ `L = K̂^{ur}`)の場合のみ**である。
原典の (i)(`𝒪/𝔭^m ≅ µ_{f,m}`)・(ii) の残り(`N_{L^m_f/L}(−α) = π^{ϕ^{m−1}}`、`α` が素元、
完全分岐)は本ファイルには**含まれない**。次数の部分(ii)は #1
(`LubinTateCompletionDegree.lean`)が、`L = K` の場合の (iii) は
`LubinTateReciprocityIsomorphism.lean` の `galoisReciprocityEquiv` が既に持っている。

## 何を示したか

**`Gal(K̂^m_f/K̂^{ur}) ≃* (𝒪_K/π^n𝒪_K)^×`**(`galoisCompletionReciprocityEquiv` /
`nonempty_galoisCompletionReciprocityEquiv`)。

## 原典との差 —— 「原典の 2 段目の同型は経由しない」

原典 (iii) は `Gal(L^m_f/L) ≅ Aut_𝒪(µ_{f,m}) ≅ (𝒪/𝔭^m)^×` と**2 段**で書く。
中央の `Aut_𝒪(µ_{f,m})`(階数 1 の自由 `𝒪/𝔭^m`-加群の自己同型群)は
(i) を前提にしており、`L = K̂^{ur}` の場合にそれを立て直すのは高い。

★★本ファイルは代わりに**基底変換の 1 段**で済ませる:

> **`Gal(K̂^{ur}(α)/K̂^{ur}) ≃* Gal(K(α)/K)`**(`galoisCompletionEquivBase`)

を作り、右辺に木の `galoisReciprocityEquiv`(`L = K` の場合の (iii)、425 行)を
**そのまま合成する**。すなわち `L = K` 版は**丸ごと流用できた**
(#1 が「次数は書き直しになった」と報告したのとは逆である)。

★★**なぜ 1 段で済むか**: `α ∈ Λ^×_{f,n}` は `𝒪_K` 係数の Eisenstein 多項式 `ψ_n` の根で、
`ψ_n` は `𝒪_{K̂^{ur}}` 上でも Eisenstein(#1 の `isEisensteinAt_map_baseIntHom_…`)だから
`K̂^{ur}` 上でも既約であり、`minpoly_{K̂^{ur}}(ι α) = ψ_n = minpoly_K(α)`。
つまり**底を `K` から `K̂^{ur}` に取り替えても最小多項式が変わらない**。
これが「不分岐拡大で底変換しても Galois 群が変わらない」の実体である。

## 何を足したか

### 抽象核(分岐・付値・Lubin-Tate の語彙が 1 つも出てこない)

| 宣言 | 内容 |
|---|---|
| `minpoly_eq_map_of_isEisensteinAt` | `A` 整閉整域・`F = Frac A`・`g ∈ A[X]` がモニック Eisenstein、`α ∈ M ⊇ F` が根なら `minpoly F α = g.map (algebraMap A F)`。#1 の `finrank_adjoin_of_isEisensteinAt` の強化版(次数だけでなく最小多項式そのもの) |
| `exists_root_preimage_of_splits` | `q ∈ L[X]` が `L` で分解するなら、体の準同型 `φ : L →+* M` について `q.map φ` の `M` での根は `φ` の像に限る |
| `map_mem_adjoin_simple_of_mem_adjoin_simple` | `φ : L →ₐ[F] M`、`z ∈ F⟮x⟯` なら `φ z ∈ F'⟮φ x⟯`(`F ⊆ F' ⊆ M` は任意の底) |
| ★★`coe_algEquiv_eq_of_adjoin_simple` | ★**基底変換の心臓**。`τ ∈ Aut(F'⟮φ x⟯/F')` と `σ : F⟮x⟯ →ₐ[F] M` が `x` で一致するなら `F⟮x⟯` **全体**で一致する。`IntermediateField.algHom_ext_of_eq_adjoin` に 2 本の `F`-代数射を渡すだけ |
| `algEquiv_adjoin_simple_ext` | `F⟮x⟯` の自己同型は `x` での値で決まる |
| `exists_algEquiv_adjoin_of_minpoly_eq` | `minpoly F α = minpoly F β` かつ `F⟮β⟯ = F⟮α⟯` なら `α ↦ β` を実現する `F⟮α⟯` の自己同型がある(`PowerBasis.equivOfMinpoly`) |
| `aeval_coe_eq_zero_of_algEquiv` | 中間体の自己同型は係数体の多項式の根を根へ写す |

### 具体層

| 宣言 | 内容 |
|---|---|
| `psiCompletion` | `ψ_n` を `K̂^{ur}` 係数で見たもの |
| `minpoly_closureCompletionCoe` | ★**`minpoly_{K̂^{ur}}(ι α) = ψ_n`**(`α ∈ Λ^×_{f,n}`) |
| `exists_torsion_of_aeval_psiCompletion` | ★**`ψ_n` の `ℂ_K` での根は `ι(Λ^×_{f,n})` に限る**(`K^{al}` で分解するから) |
| `closureCompletionCoe_mem_adjoin` | `z ∈ K⟮α⟯` なら `ι z ∈ K̂^{ur}⟮ι α⟯` |
| `algEquivBase` | ★**制限写像 `Gal(K̂^{ur}(α)/K̂^{ur}) → Gal(K(α)/K)`** |
| `coe_algEquivBase` | `τ` は `ι(K(α))` の上で `algEquivBase τ` として働く |
| `algEquivBase_mul` / `_injective` / `_surjective` | 群準同型・単射・全射 |
| `galoisCompletionEquivBase` | ★★**`Gal(K̂^{ur}(α)/K̂^{ur}) ≃* Gal(K(α)/K)`** |
| `galoisCompletionReciprocityEquiv` | ★★★**`Gal(K̂^m_f/K̂^{ur}) ≃* (𝒪_K/π^n)^×`** |
| `nonempty_galoisCompletionReciprocityEquiv` | 同上、原始点の選択を隠した形 |

## 退化の自己検査

* ★★**`α ∈ µ^×_{f,m}`(位数がちょうど `π^n`)を落とすと偽**——低い層の点 `α ∈ Λ_{n−1}` では
  `ψ_n` の根でなく、`minpoly_{K̂^{ur}}(ι α) = ψ_n` が成り立たない(次数が下がる)。
  本ファイルはどの宣言でも `hxψ : x ∈ iteratedLubinTatePsiTorsionPoints … n hn` を要求する。
* ★**`m ≥ 1`(`hn : 1 ≤ n`)を落とすと `(𝒪/𝔭^m)^×` が自明**になり、`ψ_n` も定義されない。
* ★★**`K̂^{ur}` の完備性を落とすと `𝒪_{K̂^{ur}}` が DVR でなくなる**——
  `minpoly_closureCompletionCoe` は `IsIntegrallyClosed 𝒪_{K̂^{ur}}`(DVR から)を
  Gauss の補題に使う。
* ★★**`𝔪_{𝒪_{K̂^{ur}}} = (π)`(`K̂^{ur}/K` が不分岐)を落とすと偽**——`ψ_n` が
  `𝒪_{K̂^{ur}}` 上で Eisenstein でなくなり、`K̂^{ur}` 上で既約とは限らない。
  そのとき `Gal(K̂^m_f/K̂^{ur})` は `(𝒪/𝔭^m)^×` の真部分群になりうる。
* ★★**剰余体の有限性は `𝒪_K` 側にしか使っていない**——`𝒪_{K̂^{ur}}` の剰余体は代数閉
  (＝無限)なので `[Fintype (ResidueField 𝒪_{K̂^{ur}})]` を要求する在庫は使えないが、
  本ファイルはそれを 1 つも使っていない(#1 と同じ。決定 D28 の予測どおり)。

## 逸脱(記録)

1. **`L = K̂^{ur}` に固定した。** 原典 §4.1 は「a complete unramified extension `L` of `K`」と
   一般に書くが、本ファイルは `L = K̂^{ur}` の場合だけを扱う(Proposition 4.7(ii) 以降で
   実際に使われるのはこの場合)。#1 と同じ逸脱である。
2. **原典の中間項 `Aut_𝒪(µ_{f,m})` を経由していない**(上の「原典との差」)。
   結論(両端の同型)は同じで、消費側(Prop 6.14)が使うのは両端だけである。
3. **`K̂^m_f` の定義は #9(`LubinTateTowerFIndependent.lean`)のもの**を使う。
4. `galoisCompletionReciprocityEquiv` は原始点 `α` の選択に依存する。原典の `ρ_{f,m}` は
   `α` に依らない(`Aut_𝒪(µ_{f,m})` を経由するから)が、本ファイルは `α` を固定して作る。
   「選択を隠した」形は `nonempty_galoisCompletionReciprocityEquiv`。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NNReal Valued Classical

/-! ## 0. 抽象核

分岐・付値・Lubin-Tate の語彙が 1 つも出てこない部分。原典の設定に依らない。 -/

section AbstractCore

open Polynomial in
/-- ★★★★★★★★★★★★★★★★**Eisenstein 多項式は根の最小多項式そのもの**。

`A` を整閉整域、`F = Frac A`、`M ⊇ F` を体とする。`g ∈ A[X]` がモニックで
素イデアル `𝓟` について Eisenstein なら、`M` の中の `g` の根 `α` について
`minpoly F α = g.map (algebraMap A F)`。

★#1(`LubinTateCompletionDegree.lean`)の `finrank_adjoin_of_isEisensteinAt` は
この等式から次数だけを読んでいた。基底変換には**最小多項式そのもの**が要る
(「底を `K` から `K̂^{ur}` へ替えても最小多項式が変わらない」が本ノードの要)。 -/
theorem minpoly_eq_map_of_isEisensteinAt
    {A F M : Type*} [CommRing A] [IsDomain A] [IsIntegrallyClosed A]
    [Field F] [Algebra A F] [IsFractionRing A F] [Field M] [Algebra F M]
    {𝓟 : Ideal A} (h𝓟 : 𝓟.IsPrime)
    {g : Polynomial A} (hmonic : g.Monic) (heis : g.IsEisensteinAt 𝓟)
    (hdegpos : 0 < g.natDegree)
    (α : M) (hroot : Polynomial.aeval α (g.map (algebraMap A F)) = 0) :
    minpoly F α = g.map (algebraMap A F) := by
  have hirrA : Irreducible g := heis.irreducible h𝓟 hmonic.isPrimitive hdegpos
  have hmonicF : (g.map (algebraMap A F)).Monic := hmonic.map _
  have hirrF : Irreducible (g.map (algebraMap A F)) :=
    (hmonic.irreducible_iff_irreducible_map_fraction_map (K := F)).mp hirrA
  exact (minpoly.eq_of_irreducible_of_monic hirrF hroot hmonicF).symm

/-- **分解する多項式の根は、体の拡大を大きくしても増えない**。

`q ∈ L[X]` が `L` の中で 1 次式の積に分解するなら、体の準同型 `φ : L →+* M` について
`q.map φ` の `M` での根はすべて `φ` の像である。

★これが「`ℂ_K` の中で `ψ_n` の根を探しても `K^{al}` の中の根しか出てこない」の実体。 -/
theorem exists_root_preimage_of_splits {L M : Type*} [Field L] [Field M]
    (φ : L →+* M) {q : Polynomial L} (hq0 : q ≠ 0) (hsplit : q.Splits)
    {β : M} (hβ : Polynomial.eval β (q.map φ) = 0) :
    ∃ y ∈ q.roots, φ y = β := by
  have hmap0 : q.map φ ≠ 0 := (Polynomial.map_ne_zero_iff φ.injective).mpr hq0
  have hmem : β ∈ (q.map φ).roots := by
    rw [Polynomial.mem_roots hmap0]
    exact hβ
  rw [hsplit.roots_map φ, Multiset.mem_map] at hmem
  obtain ⟨y, hy, hyb⟩ := hmem
  exact ⟨y, hy, hyb⟩

/-- **1 元生成は底を大きくしても保たれる**(元のレベル)。

`φ : L →ₐ[F] M`、`z ∈ F⟮x⟯` なら、`F` 上の任意の底 `F'` について `φ z ∈ F'⟮φ x⟯`。

★#1 の `adjoin_image_eq_adjoin_simple` の「元ごと」版。基底変換の制限写像を
作るのに、集合の等式ではなく元の所属が要る。 -/
theorem map_mem_adjoin_simple_of_mem_adjoin_simple
    {F L M F' : Type*} [Field F] [Field L] [Field M] [Field F']
    [Algebra F L] [Algebra F M] [Algebra F' M] [Algebra F F'] [IsScalarTower F F' M]
    (φ : L →ₐ[F] M) (x : L) {z : L} (hz : z ∈ IntermediateField.adjoin F ({x} : Set L)) :
    φ z ∈ IntermediateField.adjoin F' ({φ x} : Set M) := by
  have h1 : φ z ∈ (IntermediateField.adjoin F ({x} : Set L)).map φ := ⟨z, hz, rfl⟩
  rw [IntermediateField.adjoin_map, Set.image_singleton] at h1
  have h2 : IntermediateField.adjoin F ({φ x} : Set M)
      ≤ (IntermediateField.adjoin F' ({φ x} : Set M)).restrictScalars F :=
    IntermediateField.adjoin_le_iff.mpr (by
      rintro _ rfl
      exact IntermediateField.subset_adjoin F' _ rfl)
  exact h2 h1

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**基底変換の心臓**。

`F ⊆ F' ⊆ M`、`φ : L →ₐ[F] M`、`x ∈ L` とする。`F'⟮φ x⟯` の `F'`-自己同型 `τ` と
`F`-代数射 `σ : F⟮x⟯ →ₐ[F] M` が **`x` の 1 点で一致する**なら、`F⟮x⟯` **全体**で一致する:
`τ (φ z) = σ z` (`∀ z ∈ F⟮x⟯`)。

段取り: `z ↦ τ ⟨φ z⟩` は `F⟮x⟯ →ₐ[F] M` の合成(`φ` の像が `F'⟮φ x⟯` に入ることは
`map_mem_adjoin_simple_of_mem_adjoin_simple`、そこへ `AlgHom.codRestrict`)なので、
`IntermediateField.algHom_ext_of_eq_adjoin` が 2 本を 1 点で比べてくれる。

★★**分岐も付値も Lubin-Tate も出てこない。** 「不分岐拡大で底変換しても Galois 群が
変わらない」という主張の**群論的・体論的な核**がこれである。 -/
theorem coe_algEquiv_eq_of_adjoin_simple
    {F F' L M : Type*} [Field F] [Field F'] [Field L] [Field M]
    [Algebra F L] [Algebra F M] [Algebra F' M] [Algebra F F'] [IsScalarTower F F' M]
    (φ : L →ₐ[F] M) (x : L)
    (σ : ↥(IntermediateField.adjoin F ({x} : Set L)) →ₐ[F] M)
    (τ : ↥(IntermediateField.adjoin F' ({φ x} : Set M)) ≃ₐ[F']
         ↥(IntermediateField.adjoin F' ({φ x} : Set M)))
    (hsx : (τ ⟨φ x, IntermediateField.mem_adjoin_simple_self F' (φ x)⟩ : M)
      = σ ⟨x, IntermediateField.mem_adjoin_simple_self F x⟩)
    (z : L) (hz : z ∈ IntermediateField.adjoin F ({x} : Set L)) :
    (τ ⟨φ z, map_mem_adjoin_simple_of_mem_adjoin_simple φ x hz⟩ : M) = σ ⟨z, hz⟩ := by
  have hkey :
      ((IntermediateField.adjoin F' ({φ x} : Set M)).val.restrictScalars F).comp
        ((((τ : _ ≃ₐ[F'] _).restrictScalars F).toAlgHom).comp
          (((φ.comp (IntermediateField.adjoin F ({x} : Set L)).val).codRestrict
            ((IntermediateField.adjoin F' ({φ x} : Set M)).toSubalgebra.restrictScalars F)
            (fun w => map_mem_adjoin_simple_of_mem_adjoin_simple φ x w.2))))
      = σ := by
    refine IntermediateField.algHom_ext_of_eq_adjoin (F := F)
      (S := IntermediateField.adjoin F ({x} : Set L)) (s := ({x} : Set L)) rfl ?_
    intro y hy
    simp only [Set.mem_singleton_iff] at hy
    subst hy
    exact hsx
  exact congrArg (fun g => g ⟨z, hz⟩) hkey

/-- **`F⟮x⟯` の `F`-自己同型は `x` での値で決まる**。 -/
theorem algEquiv_adjoin_simple_ext {F M : Type*} [Field F] [Field M] [Algebra F M] (x : M)
    (ρ ρ' : ↥(IntermediateField.adjoin F ({x} : Set M)) ≃ₐ[F]
            ↥(IntermediateField.adjoin F ({x} : Set M)))
    (h : (ρ ⟨x, IntermediateField.mem_adjoin_simple_self F x⟩ : M)
       = (ρ' ⟨x, IntermediateField.mem_adjoin_simple_self F x⟩ : M)) : ρ = ρ' := by
  apply AlgEquiv.coe_algHom_injective
  refine IntermediateField.algHom_ext_of_eq_adjoin (F := F)
    (S := IntermediateField.adjoin F ({x} : Set M)) (s := ({x} : Set M)) rfl ?_
  intro y hy
  simp only [Set.mem_singleton_iff] at hy
  subst hy
  exact Subtype.ext h

/-- **最小多項式が同じ 2 元は、生成する体が同じなら自己同型で入れ替わる**。

`minpoly F α = minpoly F β` かつ `F⟮β⟯ = F⟮α⟯` なら、`α ↦ β` を実現する
`F⟮α⟯` の `F`-自己同型がある(`PowerBasis.equivOfMinpoly` + `IntermediateField.equivOfEq`)。 -/
theorem exists_algEquiv_adjoin_of_minpoly_eq {F M : Type*} [Field F] [Field M] [Algebra F M]
    {α β : M} (hα : IsIntegral F α) (hβ : IsIntegral F β)
    (hmin : minpoly F α = minpoly F β)
    (hFeq : IntermediateField.adjoin F ({β} : Set M) = IntermediateField.adjoin F ({α} : Set M)) :
    ∃ τ : ↥(IntermediateField.adjoin F ({α} : Set M)) ≃ₐ[F]
          ↥(IntermediateField.adjoin F ({α} : Set M)),
      (τ ⟨α, IntermediateField.mem_adjoin_simple_self F α⟩ : M) = β := by
  have hgen : minpoly F (IntermediateField.adjoin.powerBasis hα).gen
      = minpoly F (IntermediateField.adjoin.powerBasis hβ).gen := by
    rw [IntermediateField.adjoin.powerBasis_gen, IntermediateField.adjoin.powerBasis_gen,
      IntermediateField.minpoly_gen, IntermediateField.minpoly_gen, hmin]
  refine ⟨((IntermediateField.adjoin.powerBasis hα).equivOfMinpoly
    (IntermediateField.adjoin.powerBasis hβ) hgen).trans (IntermediateField.equivOfEq hFeq), ?_⟩
  show ((IntermediateField.equivOfEq hFeq)
    (((IntermediateField.adjoin.powerBasis hα).equivOfMinpoly
      (IntermediateField.adjoin.powerBasis hβ) hgen)
        ⟨α, IntermediateField.mem_adjoin_simple_self F α⟩) : M) = β
  have hg : (⟨α, IntermediateField.mem_adjoin_simple_self F α⟩ :
      ↥(IntermediateField.adjoin F ({α} : Set M)))
      = (IntermediateField.adjoin.powerBasis hα).gen := by
    rw [IntermediateField.adjoin.powerBasis_gen]; rfl
  rw [hg, PowerBasis.equivOfMinpoly_gen, IntermediateField.adjoin.powerBasis_gen]
  rfl

/-- **中間体の自己同型は、係数体の多項式の根を根へ写す**。 -/
theorem aeval_coe_eq_zero_of_algEquiv {F M : Type*} [Field F] [Field M] [Algebra F M]
    {S : IntermediateField F M} (τ : ↥S ≃ₐ[F] ↥S) (q : Polynomial F) (w : ↥S)
    (hw : Polynomial.aeval ((w : M)) q = 0) :
    Polynomial.aeval (((τ w : ↥S) : M)) q = 0 := by
  have key : ∀ v : ↥S, Polynomial.aeval ((v : M)) q = ((Polynomial.aeval v q : ↥S) : M) :=
    fun v => Polynomial.aeval_algHom_apply S.val v q
  have h1 : Polynomial.aeval w q = 0 := by
    have h2 := key w
    rw [hw] at h2
    exact Subtype.ext h2.symm
  have h2 : Polynomial.aeval (τ w) q = τ (Polynomial.aeval w q) :=
    Polynomial.aeval_algHom_apply τ w q
  rw [key (τ w), h2, h1, map_zero]
  rfl

end AbstractCore

variable {p : ℕ} [Fact p.Prime]

section Concrete

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

/-! ## 1. `ψ_n` を `K̂^{ur}` 係数で見る -/

/-- `ψ_n` を `𝒪_K → 𝒪_{K̂^{ur}} → K̂^{ur}` で写した多項式。 -/
noncomputable def psiCompletion (n : ℕ) (hn : 1 ≤ n) : Polynomial (unramifiedCompletion K) :=
  ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map (baseIntHom K)).map
    (algebraMap ↥(unramifiedCompletionInt K) (unramifiedCompletion K))

/-- ★★★★★★★★★★★★★★★★★★★★**`minpoly_{K̂^{ur}}(ι α) = ψ_n`**(`α ∈ Λ^×_{f,n}`)。

★これが本ノードの要である。`ψ_n` は `𝒪_{K̂^{ur}}` 上でも Eisenstein
(#1 の `isEisensteinAt_map_baseIntHom_iteratedLubinTatePsi`)なので `K̂^{ur}` 上でも既約で、
**底を `K` から `K̂^{ur}` へ替えても最小多項式が変わらない**。 -/
theorem minpoly_closureCompletionCoe (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    minpoly (unramifiedCompletion K) (closureCompletionCoe K x)
      = psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn := by
  haveI := isDiscreteValuationRing_unramifiedCompletionInt K
  refine minpoly_eq_map_of_isEisensteinAt
    (IsLocalRing.maximalIdeal.isMaximal ↥(unramifiedCompletionInt K)).isPrime
    ((isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic.map _)
    (isEisensteinAt_map_baseIntHom_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn)
    ?_ _
    (aeval_closureCompletionCoe_map_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ)
  rw [(isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic.natDegree_map,
    natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn]
  exact sub_pow_pred_pos K hq n hn

/-- `ψ_n` を `K̂^{ur}` 経由で `ℂ_K` へ写す道は、`K^{al}` 経由で写す道に等しい。 -/
theorem map_psiCompletion (n : ℕ) (hn : 1 ≤ n) :
    (psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn).map
        (algebraMap (unramifiedCompletion K) (closureCompletion K))
      = ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map
          (algebraMap 𝒪[K.carrier] K.closure)).map (closureCompletionCoe K) := by
  rw [psiCompletion, Polynomial.map_map, Polynomial.map_map, Polynomial.map_map,
    ← algebraMap_comp_baseIntHom_eq K]

/-- ★★★★★★★★★★★★★★★★**`ψ_n` の `ℂ_K` での根は `ι(Λ^×_{f,n})` に限る**。

`ψ_n` は `K^{al}` の中で既に 1 次式の積に分解している(`K^{al}` は代数閉)ので、
`ℂ_K` へ写しても根は増えない(抽象核 `exists_root_preimage_of_splits`)。

★★これが無いと「`τ(ι α)` は `ι` の像に入っている」が言えず、制限写像が作れない。 -/
theorem exists_torsion_of_aeval_psiCompletion (n : ℕ) (hn : 1 ≤ n)
    (β : closureCompletion K)
    (hβ : Polynomial.aeval β (psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn) = 0) :
    ∃ y ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn,
      closureCompletionCoe K y = β := by
  have hmonic := (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic
  have hq0 : (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map
      (algebraMap 𝒪[K.carrier] K.closure) ≠ 0 := (hmonic.map _).ne_zero
  have hsplit : ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map
      (algebraMap 𝒪[K.carrier] K.closure)).Splits := IsAlgClosed.splits _
  have hβ' : Polynomial.eval β
      (((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map
        (algebraMap 𝒪[K.carrier] K.closure)).map (closureCompletionCoe K)) = 0 := by
    rw [← map_psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn, ← Polynomial.eval₂_eq_eval_map,
      ← Polynomial.aeval_def]
    exact hβ
  obtain ⟨y, hy, hyβ⟩ := exists_root_preimage_of_splits (closureCompletionCoe K) hq0 hsplit hβ'
  exact ⟨y, Multiset.mem_toFinset.mpr hy, hyβ⟩

/-! ## 2. 制限写像 `Gal(K̂^{ur}(α)/K̂^{ur}) → Gal(K(α)/K)` -/

omit [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
  [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))] in
/-- `z ∈ K⟮α⟯` なら `ι z ∈ K̂^{ur}⟮ι α⟯`。 -/
theorem closureCompletionCoe_mem_adjoin (x : K.closure) {z : K.closure}
    (hz : z ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    closureCompletionCoe K z ∈ IntermediateField.adjoin (unramifiedCompletion K)
      ({closureCompletionCoe K x} : Set (closureCompletion K)) :=
  map_mem_adjoin_simple_of_mem_adjoin_simple (closureCompletionAlgHom K) x hz

/-- ★★★★★★★★★★★★★★★★★★★★**制限写像の存在と一意性**。

`τ ∈ Gal(K̂^{ur}(α)/K̂^{ur})` に対し、`ι(ρ α) = τ(ι α)` を満たす
`ρ ∈ Gal(K(α)/K)` がただ 1 つ存在する。

* 存在: `τ(ι α)` は `ψ_n` の `ℂ_K` での根(抽象核 `aeval_coe_eq_zero_of_algEquiv`)なので
  `exists_torsion_of_aeval_psiCompletion` で `ι y`(`y ∈ Λ^×_{f,n}`)の形に書け、
  木の `exists_algEquiv_of_mem_iteratedLubinTatePsiTorsionPoints` が `σ α = y` なる
  `σ ∈ Gal(K^{al}/K)` を与える。その `K⟮α⟯` への制限が `ρ`。
* 一意性: `K⟮α⟯` の自己同型は `α` での値で決まる(抽象核 `algEquiv_adjoin_simple_ext`)。 -/
theorem existsUnique_algEquivBase (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (τ : ↥(IntermediateField.adjoin (unramifiedCompletion K)
            ({closureCompletionCoe K x} : Set (closureCompletion K)))
          ≃ₐ[unramifiedCompletion K]
         ↥(IntermediateField.adjoin (unramifiedCompletion K)
            ({closureCompletionCoe K x} : Set (closureCompletion K)))) :
    ∃! ρ : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ≃ₐ[K.carrier]
             ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)),
      closureCompletionCoe K ((ρ ⟨x, hmem⟩ : ↥(IntermediateField.adjoin K.carrier
          ({x} : Set K.closure))) : K.closure)
        = (τ ⟨closureCompletionCoe K x,
              IntermediateField.mem_adjoin_simple_self (unramifiedCompletion K) _⟩ :
            closureCompletion K) := by
  have hroot : Polynomial.aeval
      ((τ ⟨closureCompletionCoe K x,
        IntermediateField.mem_adjoin_simple_self (unramifiedCompletion K) _⟩ :
          closureCompletion K)) (psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn) = 0 :=
    aeval_coe_eq_zero_of_algEquiv τ _ _
      (aeval_closureCompletionCoe_map_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ)
  obtain ⟨y, hy, hyeq⟩ :=
    exists_torsion_of_aeval_psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn _ hroot
  obtain ⟨s, hs⟩ := exists_algEquiv_of_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0
    f hf0 hf1 hf n hn x y hxψ hy
  have hmain : closureCompletionCoe K
      ((algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem s ⟨x, hmem⟩ :
        ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure)
      = (τ ⟨closureCompletionCoe K x,
            IntermediateField.mem_adjoin_simple_self (unramifiedCompletion K) _⟩ :
          closureCompletion K) := by
    rw [coe_algEquivRestrictSelf]
    show closureCompletionCoe K (s x) = _
    rw [hs]
    exact hyeq
  refine ⟨algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem s, hmain, ?_⟩
  · intro ρ hρ
    refine algEquiv_adjoin_simple_ext x ρ _ ?_
    apply closureCompletionCoe_injective K
    rw [hρ, coe_algEquivRestrictSelf]
    show _ = closureCompletionCoe K (s x)
    rw [hs]
    exact hyeq.symm

/-- **制限写像 `Gal(K̂^{ur}(α)/K̂^{ur}) → Gal(K(α)/K)`**。 -/
noncomputable def algEquivBase (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (τ : ↥(IntermediateField.adjoin (unramifiedCompletion K)
            ({closureCompletionCoe K x} : Set (closureCompletion K)))
          ≃ₐ[unramifiedCompletion K]
         ↥(IntermediateField.adjoin (unramifiedCompletion K)
            ({closureCompletionCoe K x} : Set (closureCompletion K)))) :
    ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ≃ₐ[K.carrier]
      ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
  (existsUnique_algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ).choose

theorem algEquivBase_spec (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (τ : ↥(IntermediateField.adjoin (unramifiedCompletion K)
            ({closureCompletionCoe K x} : Set (closureCompletion K)))
          ≃ₐ[unramifiedCompletion K]
         ↥(IntermediateField.adjoin (unramifiedCompletion K)
            ({closureCompletionCoe K x} : Set (closureCompletion K)))) :
    closureCompletionCoe K
        ((algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ ⟨x, hmem⟩ :
          ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure)
      = (τ ⟨closureCompletionCoe K x,
            IntermediateField.mem_adjoin_simple_self (unramifiedCompletion K) _⟩ :
          closureCompletion K) :=
  (existsUnique_algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ).choose_spec.1

/-- ★★★★★★★★★★★★★★**`τ` は `ι(K(α))` の上で `algEquivBase τ` として働く**。

`algEquivBase` は `α` の 1 点でしか定義していないが、抽象核
`coe_algEquiv_eq_of_adjoin_simple` によって `K⟮α⟯` 全体に伝播する。
★群準同型性(`algEquivBase_mul`)はここから出る。 -/
theorem coe_algEquivBase (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (τ : ↥(IntermediateField.adjoin (unramifiedCompletion K)
            ({closureCompletionCoe K x} : Set (closureCompletion K)))
          ≃ₐ[unramifiedCompletion K]
         ↥(IntermediateField.adjoin (unramifiedCompletion K)
            ({closureCompletionCoe K x} : Set (closureCompletion K))))
    (z : K.closure) (hz : z ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    (τ ⟨closureCompletionCoe K z, closureCompletionCoe_mem_adjoin K x hz⟩ :
        closureCompletion K)
      = closureCompletionCoe K
          ((algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ ⟨z, hz⟩ :
            ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure) :=
  coe_algEquiv_eq_of_adjoin_simple (F' := unramifiedCompletion K) (closureCompletionAlgHom K) x
    (((closureCompletionAlgHom K).comp
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).val).comp
      (algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ).toAlgHom)
    τ (algEquivBase_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ).symm z hz

/-- **`algEquivBase` は群準同型**。 -/
theorem algEquivBase_mul (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (τ τ' : ↥(IntermediateField.adjoin (unramifiedCompletion K)
            ({closureCompletionCoe K x} : Set (closureCompletion K)))
          ≃ₐ[unramifiedCompletion K]
         ↥(IntermediateField.adjoin (unramifiedCompletion K)
            ({closureCompletionCoe K x} : Set (closureCompletion K)))) :
    algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem (τ * τ')
      = algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ *
        algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ' := by
  refine (existsUnique_algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem (τ * τ')).unique
    (algEquivBase_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem (τ * τ')) ?_
  have h1 : (τ' ⟨closureCompletionCoe K x,
        IntermediateField.mem_adjoin_simple_self (unramifiedCompletion K) _⟩)
      = ⟨closureCompletionCoe K
          ((algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ' ⟨x, hmem⟩ :
            ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure),
          closureCompletionCoe_mem_adjoin K x
            (algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ' ⟨x, hmem⟩).2⟩ :=
    Subtype.ext
      (algEquivBase_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ').symm
  rw [AlgEquiv.mul_apply, AlgEquiv.mul_apply, h1]
  exact (coe_algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ _
    (algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ' ⟨x, hmem⟩).2).symm

/-- **`algEquivBase` は単射**——`K̂^{ur}(α)` の自己同型は `ι α` での値で決まるから。 -/
theorem algEquivBase_injective (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Function.Injective
      (algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem) := by
  intro τ τ' h
  refine algEquiv_adjoin_simple_ext (closureCompletionCoe K x) τ τ' ?_
  rw [← algEquivBase_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ,
    ← algEquivBase_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ', h]

/-- ★★★★★★★★★★★★★★★★**`algEquivBase` は全射**。

`ρ ∈ Gal(K(α)/K)` に対し `y := ρ α` はまた `ψ_n` の根なので、
`minpoly_{K̂^{ur}}(ι y) = ψ_n = minpoly_{K̂^{ur}}(ι α)` であり、しかも
`K̂^{ur}(ι y) = K̂^m_f = K̂^{ur}(ι α)`(#1 の `lubinTateCompletionField_eq_adjoin_simple`)。
抽象核 `exists_algEquiv_adjoin_of_minpoly_eq` が `ι α ↦ ι y` を実現する `τ` を与える。 -/
theorem algEquivBase_surjective (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Function.Surjective
      (algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem) := by
  intro ρ
  have hy : ((ρ ⟨x, hmem⟩ : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
      K.closure) ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn :=
    algEquivL_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hmem ρ
  have hyn : ((ρ ⟨x, hmem⟩ : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
      K.closure) ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n := by
    rw [iteratedLubinTateTorsionPoints_eq_union K hq hπmax hπne0 f hf0 hf1 hf n hn,
      Finset.mem_union]
    exact Or.inr hy
  have hmin : minpoly (unramifiedCompletion K) (closureCompletionCoe K x)
      = minpoly (unramifiedCompletion K) (closureCompletionCoe K
        ((ρ ⟨x, hmem⟩ : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
          K.closure)) := by
    rw [minpoly_closureCompletionCoe K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ,
      minpoly_closureCompletionCoe K hq hπmax hπne0 f hf0 hf1 hf n hn _ hy]
  have hadj : IntermediateField.adjoin (unramifiedCompletion K)
        ({closureCompletionCoe K
          ((ρ ⟨x, hmem⟩ : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
            K.closure)} : Set (closureCompletion K))
      = IntermediateField.adjoin (unramifiedCompletion K)
        ({closureCompletionCoe K x} : Set (closureCompletion K)) :=
    (lubinTateCompletionField_eq_adjoin_simple K hq hπmax hπne0 f hf0 hf1 hf n hn _ hy hyn).symm.trans
      (lubinTateCompletionField_eq_adjoin_simple K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn)
  obtain ⟨τ, hτ⟩ := exists_algEquiv_adjoin_of_minpoly_eq
    (isIntegral_closureCompletionCoe K x) (isIntegral_closureCompletionCoe K _) hmin hadj
  refine ⟨τ, ?_⟩
  exact (existsUnique_algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ).unique
    (algEquivBase_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ) hτ.symm

/-! ## 3. 主結果 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Gal(K̂^{ur}(α)/K̂^{ur}) ≃* Gal(K(α)/K)`** —— 不分岐な底変換で Galois 群は変わらない。

★これが本ノードの主張の実体である。原典 Proposition 4.4(iii) の `L = K̂` の場合は、
これと `L = K` の場合(木の `galoisReciprocityEquiv`)の合成で出る。 -/
noncomputable def galoisCompletionEquivBase (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    (↥(IntermediateField.adjoin (unramifiedCompletion K)
        ({closureCompletionCoe K x} : Set (closureCompletion K)))
      ≃ₐ[unramifiedCompletion K]
     ↥(IntermediateField.adjoin (unramifiedCompletion K)
        ({closureCompletionCoe K x} : Set (closureCompletion K)))) ≃*
    (↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ≃ₐ[K.carrier]
     ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
  MulEquiv.ofBijective
    (MonoidHom.mk' (algEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem)
      (algEquivBase_mul K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem))
    ⟨algEquivBase_injective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem,
      algEquivBase_surjective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem⟩

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Gal(K̂^m_f/K̂^{ur}) ≃* (𝒪_K/π^n𝒪_K)^×`**
(Yoshida 2008 Proposition 4.4(iii)、`L = K̂` の場合)。

3 本の合成:
`Gal(K̂^m_f/K̂^{ur}) ≃ Gal(K̂^{ur}(α)/K̂^{ur})`(#1 の `lubinTateCompletionField_eq_adjoin_simple`)
`≃ Gal(K(α)/K)`(`galoisCompletionEquivBase`)
`≃ (𝒪_K/π^n)^×`(木の `galoisReciprocityEquiv`)。 -/
noncomputable def galoisCompletionReciprocityEquiv (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    (↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n) ≃ₐ[unramifiedCompletion K]
      ↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n)) ≃*
      (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))ˣ :=
  ((AlgEquiv.autCongr (IntermediateField.equivOfEq
      (lubinTateCompletionField_eq_adjoin_simple K hq hπmax hπne0 f hf0 hf1 hf n hn x
        hxψ hxn))).trans
    (galoisCompletionEquivBase K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem)).trans
      (galoisReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Gal(K̂^m_f/K̂^{ur}) ≃* (𝒪_K/π^n𝒪_K)^×`**(原始点の選択を隠した形)。

原典 Proposition 4.4(iii) の `L = K̂`(＝`L = K̂^{ur}`)の場合。
原始点 `α ∈ Λ^×_{f,n}` の存在は `iteratedLubinTatePsiTorsionPoints_nonempty`(木に在る)から。 -/
theorem nonempty_galoisCompletionReciprocityEquiv (n : ℕ) (hn : 1 ≤ n) :
    Nonempty ((↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n)
        ≃ₐ[unramifiedCompletion K]
      ↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n)) ≃*
      (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))ˣ) := by
  obtain ⟨x, hxψ⟩ :=
    iteratedLubinTatePsiTorsionPoints_nonempty K hq hπmax hπne0 f hf0 hf1 hf n hn
  have hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n := by
    rw [iteratedLubinTateTorsionPoints_eq_union K hq hπmax hπne0 f hf0 hf1 hf n hn,
      Finset.mem_union]
    exact Or.inr hxψ
  haveI := finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0
    f hf0 hf1 hf n x hxn
  exact ⟨galoisCompletionReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn
    (IntermediateField.mem_adjoin_simple_self K.carrier x)⟩

def nonempty_galoisCompletionReciprocityEquiv.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 7, item := "Proposition 4.4", sectionId := "prop-4-4" }

end Concrete

end ABC3.Found.PGC
