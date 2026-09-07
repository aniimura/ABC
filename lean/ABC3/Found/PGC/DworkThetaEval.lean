import ABC3.Found.PGC.ClosureCompletion
import ABC3.Found.PGC.DworkThetaStep2
import ABC3.Found.PGC.LubinTateFieldFIndependent

/-!
# `θ` が捩れ点の全単射 `Λ_{f,n} → Λ_{g,n}` を誘導する(`sorry` 無し)

経路 Λ の節点 **#7**。Λ6(`DworkThetaStep2.lean`)が作った

```
(σθ) ∘ f = g ∘ θ      (f ∈ F_π, g ∈ F_ϖ, ϖ = uπ)
```

と、M3(`ClosureCompletion.lean`)が作った `θ(λ)` の住処 `𝒪_{ℂ_K}` を繋いで、
原典 Milne CFT が Prop. 3.10 の**消費側**で行う議論

```
f_n(λ) = 0  ⟹  g_n(θ(λ)) = 0  ⟹  θ(λ) ∈ ι(Λ_{g,n})
```

を形式化し、さらに `θ` が `Λ_{f,n}` から `Λ_{g,n}` への**全単射**を誘導することを示す。

## 典拠(Milne, `Class Field Theory`, v4.03)

物理ページ 51、`PROOF (THAT K_π · K^un IS INDEPENDENT OF π)` の中:

> From Proposition 3.10 we find that
> (σθ) ∘ [π]_f = θ ∘ [u]_f ∘ [π]_f = θ ∘ [ϖ]_f = [ϖ]_g ∘ θ,
> that is, that (σθ)(f(T)) = g(θ(T)).
> Therefore, for any α ∈ K^al (recall that this is the separable algebraic closure of K),
> f_n(α) = 0 ⟹ g_n(θ(α)) = 0, and, similarly, g_n(α) = 0 ⟹ f_n(θ^{-1}(α)) = 0.
> Therefore θ defines a bijection Λ_{f,1} → Λ_{g,1}

および、同じ証明の末尾:

> The argument extends without difficulty to show that K^un[Λ_{g,n}] = K^un[Λ_{f,n}]
> for all n

★引用は `ResearchPaper/0_Source/Milne - Class Field Theory.txt`(pdftotext 出力)
から取ったが、pdftotext は**下付き添字と `Λ` を落とす**(原文の
`f_n(α) = 0 ⟹ g_n(θ(α)) = 0` は `f .˛/ D 0 ) g..˛// D 0` と出る)。
上の引用は添字と記号を**補って**読める形にしてある。

★本ファイルは 2 つめの引用の合図語 `without difficulty` を潰す:
原典が `n = 1` で書いて「同様に一般の `n` でも」と畳んだところを、
はじめから**任意の `n`** について証明している。

★`.src` は**書いていない**。典拠の Milne CFT は
`ResearchPaper/1_Structured/` に構造化されておらず、`ABC3.Meta.Source` の
必須フィールド `sectionId` を正直に埋められないため(嘘の `sectionId` は
書かない)。Λ6a′(`LubinTateFieldFIndependent.lean`)・M3
(`ClosureCompletion.lean`)と同じ扱いで、逐語引用と物理ページを
この docstring に置くに留める。

## 何を足したか

### 抽象核(分岐・付値・Galois の語彙が 1 つも出てこない)

一般の可換環 `A` 上の形式冪級数と、一般の整域上の多項式だけの主張:

| 宣言 | 内容 |
|---|---|
| `constantCoeff_iterate_map_eq_zero` / `hasSubst_iterate_map` | 係数写像を `n` 回当てても定数項 `0` は保たれる |
| `iterate_map_subst` | `ψ^n(P ∘ a) = ψ^n(P) ∘ ψ^n(a)` |
| `subst_iterate_map_intertwine` | `(σ^{m+1}θ) ∘ f = g ∘ (σ^m θ)`(σ をずらした絡み) |
| `subst_iteratedLubinTate_of_semilinear_intertwine` | ★★**#7a の核**: `(σθ)∘f = g∘θ` ⟹ `(σ^nθ)∘f_n = g_n∘θ` |
| `subst_intertwine_of_comp_inverse` | 合成逆 `θ'` は絡みを逆向きに満たす: `(σθ')∘g = f∘θ'` |
| `map_iteratedLubinTate` | 係数写像は自己合成と可換 |
| `exists_mem_roots_of_eval_map_eq_zero` | ★★**根の同定の核**: 根の個数が次数に等しいモニック多項式の像の根は、もとの根の像(**代数閉性は要らない。整域性だけ**) |

★`subst_iteratedLubinTate_of_semilinear_intertwine` は `CommRing A` と
`ψ : A →+* A` だけで成り立つ。Λ6a′ の `subst_iteratedLubinTate_of_intertwine`
(`ψ = id` の場合)の**半線形版**にあたる。

### 具体層

| 宣言 | 内容 |
|---|---|
| `map_baseIntHom_fixed` | `𝒪_K` 係数の級数は `σ` で動かない |
| `subst_iteratedLubinTate_of_dworkTheta` | ★**#7a**(抽象核への代入) |
| `norm_aeval_closureCompletionInt_le` 他 | `‖θ(λ)‖ ≤ ‖λ‖ < 1`、値もまた位相的冪零 |
| `coe_algebraMap_base_closureCompletionInt` | `𝒪_K → 𝒪_{ℂ_K}` は `𝒪_K → K^al → ℂ_K` と一致 |
| `aeval_map_baseIntHom_coe` | `𝒪_K` 係数多項式の評価はスカラー塔で不変 |
| `aeval_map_baseIntHom_iteratedLubinTate_eq_zero_iff` | `[π^n]_f(z) = 0 ↔ D_n(z) = 0`(`𝒪_{ℂ_K}` 版) |
| `exists_mem_torsionPoints_of_aeval_eq_zero` | `D_n(y) = 0` ⟹ `y ∈ ι(Λ_n)` |
| `exists_mem_torsionPoints_evalAt` | ★★**#7b**: `θ(Λ_{f,n}) ⊆ ι(Λ_{g,n})` |
| `exists_mem_torsionPoints_evalAt_inverse` | 原典の「and, similarly」: `θ'(Λ_{g,n}) ⊆ ι(Λ_{f,n})` |
| `coe_aeval_evalAt_comp_inverse` | `θ'(θ(λ)) = λ` |
| `exists_bijOn_torsionPoints_of_dworkTheta` | ★★**#7c**: `θ` の誘導する `Λ_{f,n} ≃ Λ_{g,n}` |
| `exists_bijOn_iteratedLubinTateTorsionPoints` | ★★★上の `θ` を Λ6 から取ってきた**完成形** |

## 代数閉性は使っていない

`ℂ_K`(= `K^al` の完備化)が代数閉であることは**一度も要求していない**
(M3 も作っていない)。根の同定は

* `Λ_n` が `𝒪_K` 係数モニック多項式 `D_n` の根集合であること、
* `D_n` の `K^al` での根の個数がちょうど `deg D_n = q^n` であること
  (`card_roots_iteratedLubinTateDistinguished_map`)、
* `ℂ_K` が整域であること

の 3 つだけで済む(`exists_mem_roots_of_eval_map_eq_zero`)。すなわち
`∏_{r ∈ Λ_{g,n}} (θλ − ι r) = 0` から**整域性だけ**で所属が出る。

## 退化の自己検査

* `‖λ‖ < 1` を落とすと `PowerSeries.HasEval` が壊れて `θ(λ)` が定義できない
  (`hasEval_coe_of_norm_lt_one` は `<` が必須。`λ = 1` で `λⁿ ↛ 0`)。
  `Λ_n` の元については `norm_lt_one_of_mem_iteratedLubinTateTorsionPoints`
  が保証する。
* `θ` の定数項が `0` であることを落とすと `f_n(λ) = 0 ⟹ (σ^nθ)(f_n(λ)) = 0`
  の段が壊れる(`aeval_zero_eq_zero`)。
* ★`f` と `g` が**同じ素元の塔に属する**ことを落とすと偽。本ファイルは
  `f ∈ F_π`・`g ∈ F_ϖ` を別々の `hπmax`・`hϖmax` で受け取っており、
  完成形 `exists_bijOn_iteratedLubinTateTorsionPoints` では
  `ϖ = uπ`(`u` は単数)、すなわち `span {ϖ} = span {π}` に固定している。
  無関係な 2 つの素元では `θ` の存在(Λ6)がそもそも言えない。
* θ′(合成逆)を落とすと単射性が出ない(`coe_aeval_evalAt_comp_inverse` が
  `subst θ θ' = X` を使う)。`θ` が単射でなければ `|Λ_f| = |Λ_g|` から
  全射も出ない。

## 逸脱(記録)

* **全単射の出し方**。原典は 2 つの包含(`θ` と `θ^{-1}`)から全単射を
  結論するが、本ファイルは「`θ` による写像が単射(`θ'` を評価して復元)」と
  `|Λ_{f,n}| = |Λ_{g,n}| = q^n`(`card_iteratedLubinTateTorsionPoints`)から
  全射を出している。原典の逆向きの包含も
  `exists_mem_torsionPoints_evalAt_inverse` として別途証明してあるので、
  原典の議論そのものも再現できる。
* **`θ` の値の住処**。原典は `θ(α)` がどの環に住むかを述べない。本ファイルは
  M3 の `𝒪_{ℂ_K}` を使い、`Λ_n` の元との一致は `ι_al : K^al ↪ ℂ_K` の像として
  述べる(`closureCompletionCoe K r = θ(λ)`)。`ι_al` は単射なので
  情報は落ちていない。
* **`Λ_n` の定義**。原典の `Λ_{f,n}` は `[π^n]_f` の零点集合だが、本ファイルは
  在庫の `iteratedLubinTateTorsionPoints`(Weierstrass 分解
  `[π^n]_f = D_n · U_n` の `D_n` の根集合)を使う。両者が一致することは
  `aeval_map_baseIntHom_iteratedLubinTate_eq_zero_iff`(単数 `U_n` を約す)で
  本ファイル内で確認している。
* **`e` の作り方**。全単射を与える写像 `e : K^al → K^al` は選択公理
  (`choose`)で作っており、`Λ_{f,n}` の外では意味を持たない(`0` を返す)。
  `Set.BijOn` は `↑Λ_{f,n}` 上でのみ主張する。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NNReal Valued Classical

/-! ### このファイル限りのインスタンス

`𝒪_{ℂ_K}` の線形位相と連続スカラー倍を**このファイル限りの**インスタンスにする。
`ClosureCompletion.lean` はこれらを `theorem` として持っており、`closureCompletionEval`
の中では `haveI` で借りている。本ファイルは `PowerSeries.aeval` を**statement に**
出すので、`haveI`(証明の中)では間に合わない。 -/

attribute [local instance] isLinearTopology_closureCompletionInt continuousSMul_closureCompletionInt

/-! ## 0. 抽象核 —— 一般の可換環上の形式冪級数だけの主張

`PowerSeries.subst a P` は mathlib の約束で `P(a)`、すなわち `P ∘ a` である。
以下ではこの向きを固定して使う。`ψ : A →+* A` は原典の `σ`(Frobenius)にあたるが、
ここでは**ただの環準同型**でよい。 -/

section AbstractCore

variable {A B : Type*} [CommRing A] [CommRing B]

/-- 係数写像を `n` 回当てても定数項 `0` は保たれる。 -/
theorem constantCoeff_iterate_map_eq_zero (ψ : A →+* A) {a : PowerSeries A}
    (ha0 : PowerSeries.constantCoeff a = 0) (n : ℕ) :
    PowerSeries.constantCoeff ((PowerSeries.map ψ)^[n] a) = 0 := by
  induction n with
  | zero => simpa using ha0
  | succ m ih => rw [Function.iterate_succ_apply']; exact constantCoeff_map_eq_zero ψ ih

/-- 上を `PowerSeries.HasSubst` の形で言い直したもの。 -/
theorem hasSubst_iterate_map (ψ : A →+* A) {a : PowerSeries A}
    (ha0 : PowerSeries.constantCoeff a = 0) (n : ℕ) :
    PowerSeries.HasSubst ((PowerSeries.map ψ)^[n] a) :=
  PowerSeries.HasSubst.of_constantCoeff_zero' (constantCoeff_iterate_map_eq_zero ψ ha0 n)

/-- 係数写像の反復は合成と可換: `ψ^n(P ∘ a) = ψ^n(P) ∘ ψ^n(a)`。 -/
theorem iterate_map_subst (ψ : A →+* A) {a F : PowerSeries A}
    (ha0 : PowerSeries.constantCoeff a = 0) (n : ℕ) :
    (PowerSeries.map ψ)^[n] (PowerSeries.subst a F)
      = PowerSeries.subst ((PowerSeries.map ψ)^[n] a) ((PowerSeries.map ψ)^[n] F) := by
  induction n with
  | zero => simp
  | succ m ih =>
    rw [Function.iterate_succ_apply', ih,
      Function.iterate_succ_apply' (f := PowerSeries.map ψ) (n := m) a,
      Function.iterate_succ_apply' (f := PowerSeries.map ψ) (n := m) F]
    exact map_subst_powerSeries ψ (hasSubst_iterate_map ψ ha0 m)

/-- **`σ` をずらした絡み**: `(σ^{m+1}θ) ∘ f = g ∘ (σ^m θ)`。

もとの絡み `(σθ) ∘ f = g ∘ θ` の両辺に `σ^m` を当て、`f`・`g` が
`σ` で動かない(`hfψ`・`hgψ`)ことを使うだけ。 -/
theorem subst_iterate_map_intertwine (ψ : A →+* A) {f g θ : PowerSeries A}
    (hf0 : PowerSeries.constantCoeff f = 0) (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hfψ : PowerSeries.map ψ f = f) (hgψ : PowerSeries.map ψ g = g)
    (hint : PowerSeries.subst f (PowerSeries.map ψ θ) = PowerSeries.subst θ g) (m : ℕ) :
    PowerSeries.subst f ((PowerSeries.map ψ)^[m + 1] θ)
      = PowerSeries.subst ((PowerSeries.map ψ)^[m] θ) g := by
  have h := congrArg (PowerSeries.map ψ)^[m] hint
  rw [iterate_map_subst ψ hf0 m, iterate_map_subst ψ hθ0 m,
    Function.iterate_fixed hfψ m, Function.iterate_fixed hgψ m] at h
  rw [Function.iterate_succ_apply]
  exact h

/-- ★★★★★★★★★★★★★★★★**抽象核(#7a)**——半線形な絡み作用素の反復:
`(σθ) ∘ f = g ∘ θ` ならば `(σ^n θ) ∘ f^{(n)} = g^{(n)} ∘ θ`。

`PowerSeries.subst h α = α(h)` の規約で、`subst f (map ψ θ) = subst θ g` は
`(σθ)(f(T)) = g(θ(T))` を意味する。証明は代入の結合律
(`subst_subst_eq`)と `subst_iterate_map_intertwine` による帰納法だけ——
`CommRing A` と環準同型 `ψ` 以外に何も要らない(局所環も付値も分岐も
Galois も出てこない)。

Λ6a′ の `subst_iteratedLubinTate_of_intertwine` は本補題の `ψ = id` の場合に
あたる(あちらは `map ψ` が消えるので `hfψ`・`hgψ` が不要)。 -/
theorem subst_iteratedLubinTate_of_semilinear_intertwine (ψ : A →+* A) {f g θ : PowerSeries A}
    (hf0 : PowerSeries.constantCoeff f = 0) (hg0 : PowerSeries.constantCoeff g = 0)
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hfψ : PowerSeries.map ψ f = f) (hgψ : PowerSeries.map ψ g = g)
    (hint : PowerSeries.subst f (PowerSeries.map ψ θ) = PowerSeries.subst θ g) (n : ℕ) :
    PowerSeries.subst (iteratedLubinTate f n) ((PowerSeries.map ψ)^[n] θ)
      = PowerSeries.subst θ (iteratedLubinTate g n) := by
  have hHSf : PowerSeries.HasSubst f := PowerSeries.HasSubst.of_constantCoeff_zero' hf0
  have hHSθ : PowerSeries.HasSubst θ := PowerSeries.HasSubst.of_constantCoeff_zero' hθ0
  induction n with
  | zero =>
    show PowerSeries.subst PowerSeries.X θ = PowerSeries.subst θ PowerSeries.X
    rw [PowerSeries.X_subst, PowerSeries.subst_X hHSθ]
  | succ m ih =>
    show PowerSeries.subst (PowerSeries.subst (iteratedLubinTate f m) f)
        ((PowerSeries.map ψ)^[m + 1] θ)
      = PowerSeries.subst θ (PowerSeries.subst (iteratedLubinTate g m) g)
    rw [← subst_subst_eq hHSf (hasSubst_iteratedLubinTate hf0 m),
      subst_iterate_map_intertwine ψ hf0 hθ0 hfψ hgψ hint m,
      subst_subst_eq (hasSubst_iterate_map ψ hθ0 m) (hasSubst_iteratedLubinTate hf0 m),
      ih, ← subst_subst_eq (hasSubst_iteratedLubinTate hg0 m) hHSθ]

/-- ★★★★★★★★**抽象核**——合成逆 `θ'` は絡みを逆向きに満たす:
`(σθ) ∘ f = g ∘ θ` かつ `θ ∘ θ' = θ' ∘ θ = X` ならば `(σθ') ∘ g = f ∘ θ'`。

原典 51 ページの「and, similarly, `g_n(α) = 0 ⟹ f_n(θ^{-1}(α)) = 0`」の
「similarly」の中身。段取りは
1. `hint` の右から `θ'` を合成して `(σθ) ∘ (f ∘ θ') = g`、
2. 左から `σθ'` を合成し、`σθ' ∘ σθ = σ(θ' ∘ θ) = σX = X` で潰す。

★`σ` は**逆にしない**(`σ^{-1}` は出てこない)。 -/
theorem subst_intertwine_of_comp_inverse (ψ : A →+* A) {f g θ θ' : PowerSeries A}
    (hf0 : PowerSeries.constantCoeff f = 0)
    (hθ0 : PowerSeries.constantCoeff θ = 0) (hθ'0 : PowerSeries.constantCoeff θ' = 0)
    (hθθ' : PowerSeries.subst θ' θ = PowerSeries.X)
    (hθ'θ : PowerSeries.subst θ θ' = PowerSeries.X)
    (hint : PowerSeries.subst f (PowerSeries.map ψ θ) = PowerSeries.subst θ g) :
    PowerSeries.subst g (PowerSeries.map ψ θ') = PowerSeries.subst θ' f := by
  have hHSf : PowerSeries.HasSubst f := PowerSeries.HasSubst.of_constantCoeff_zero' hf0
  have hHSθ : PowerSeries.HasSubst θ := PowerSeries.HasSubst.of_constantCoeff_zero' hθ0
  have hHSθ' : PowerSeries.HasSubst θ' := PowerSeries.HasSubst.of_constantCoeff_zero' hθ'0
  have hHSmθ : PowerSeries.HasSubst (PowerSeries.map ψ θ) :=
    PowerSeries.HasSubst.of_constantCoeff_zero' (constantCoeff_map_eq_zero ψ hθ0)
  have hHSθ'f : PowerSeries.HasSubst (PowerSeries.subst θ' f) :=
    PowerSeries.HasSubst.of_constantCoeff_zero' (PowerSeries.constantCoeff_subst_eq_zero hθ'0 f hf0)
  have h1 : PowerSeries.subst (PowerSeries.subst θ' f) (PowerSeries.map ψ θ) = g := by
    have := congrArg (fun t => PowerSeries.subst θ' t) hint
    simpa [subst_subst_eq hHSf hHSθ', subst_subst_eq hHSθ hHSθ', hθθ',
      PowerSeries.X_subst] using this
  have h2 := congrArg (fun t => PowerSeries.subst t (PowerSeries.map ψ θ')) h1
  rw [← subst_subst_eq hHSmθ hHSθ'f,
    ← map_subst_powerSeries ψ hHSθ, hθ'θ, PowerSeries.map_X,
    PowerSeries.subst_X hHSθ'f] at h2
  exact h2.symm

/-- 係数写像は自己合成と可換: `φ(f^{(n)}) = (φ f)^{(n)}`。 -/
theorem map_iteratedLubinTate (φ : A →+* B) {f : PowerSeries A}
    (hf0 : PowerSeries.constantCoeff f = 0) (n : ℕ) :
    PowerSeries.map φ (iteratedLubinTate f n) = iteratedLubinTate (PowerSeries.map φ f) n := by
  induction n with
  | zero => simp [iteratedLubinTate]
  | succ m ih =>
    show PowerSeries.map φ (PowerSeries.subst (iteratedLubinTate f m) f)
      = PowerSeries.subst (iteratedLubinTate (PowerSeries.map φ f) m) (PowerSeries.map φ f)
    rw [map_subst_powerSeries φ (hasSubst_iteratedLubinTate hf0 m), ih]

/-- ★★★★★★★★★★★★**抽象核(根の同定)**——根の個数が次数に等しい
モニック多項式 `Q ∈ L[X]` と環準同型 `g : L →+* M`(`M` は整域)について、
`Q^g` の根 `y` は必ず `Q` の根の像である。

`Q = ∏_{r ∈ Q.roots} (X − r)` を `g` で写して
`∏_{r} (y − g r) = 0` とし、**整域性だけ**で 1 つの因子が `0` だと結論する。
★`L` の代数閉性も `M` の代数閉性も要らない(`hcard` が splits の代わり)。 -/
theorem exists_mem_roots_of_eval_map_eq_zero {L M : Type*} [CommRing L] [IsDomain L]
    [CommRing M] [IsDomain M] (g : L →+* M) {Q : Polynomial L} (hQ : Q.Monic)
    (hcard : Multiset.card Q.roots = Q.natDegree) {y : M} (hy : (Q.map g).eval y = 0) :
    ∃ r ∈ Q.roots, y = g r := by
  rw [← Polynomial.prod_multiset_X_sub_C_of_monic_of_roots_card_eq hQ hcard,
    Polynomial.map_multiset_prod, Multiset.map_map, Polynomial.eval_multiset_prod,
    Multiset.map_map] at hy
  simp only [Function.comp_apply, Polynomial.map_sub, Polynomial.map_X, Polynomial.map_C,
    Polynomial.eval_sub, Polynomial.eval_X, Polynomial.eval_C] at hy
  obtain ⟨r, hr, hrw⟩ := Multiset.mem_map.mp (Multiset.prod_eq_zero_iff.mp hy)
  exact ⟨r, hr, sub_eq_zero.mp hrw⟩

end AbstractCore

/-! ## 1. #7a の具体層 —— 抽象核へ代入するだけ -/

variable {p : ℕ} [Fact p.Prime]

/-- `𝒪_K` 係数の冪級数は `σ ∈ Gal(K^ur/K)` で動かない
(`unramGalCompletionInt_baseIntHom` を係数ごとに当てるだけ)。 -/
theorem map_baseIntHom_fixed (K : PAdicLocalField p) (σ : unramGal K)
    (f : PowerSeries (𝒪[K.carrier])) :
    PowerSeries.map (unramGalCompletionIntHom K σ) (PowerSeries.map (baseIntHom K) f)
      = PowerSeries.map (baseIntHom K) f := by
  refine PowerSeries.ext fun n => ?_
  simp only [PowerSeries.coeff_map]
  exact unramGalCompletionInt_baseIntHom K σ _

/-- ★★★★★★★★★★★★★★★★**#7a**——Λ6 が出す絡み
`(σθ) ∘ f = g ∘ θ` から、その `n` 回反復
`(σ^n θ) ∘ f^{(n)} = g^{(n)} ∘ θ` を得る。

抽象核 `subst_iteratedLubinTate_of_semilinear_intertwine` に
`A := 𝒪_{K̂^{ur}}`・`ψ := σ` を代入し、`f`・`g` が `σ` で動かないことを
`map_baseIntHom_fixed` で供給するだけ。★M3(`ℂ_K`)には依存しない。 -/
theorem subst_iteratedLubinTate_of_dworkTheta (K : PAdicLocalField p)
    (σ : unramGal K) {f g : PowerSeries (𝒪[K.carrier])}
    {θ : PowerSeries ↥(unramifiedCompletionInt K)}
    (hf0 : PowerSeries.constantCoeff f = 0) (hg0 : PowerSeries.constantCoeff g = 0)
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hint : PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
        (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g)) (n : ℕ) :
    PowerSeries.subst (PowerSeries.map (baseIntHom K) (iteratedLubinTate f n))
        ((PowerSeries.map (unramGalCompletionIntHom K σ))^[n] θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) (iteratedLubinTate g n)) := by
  rw [map_iteratedLubinTate _ hf0, map_iteratedLubinTate _ hg0]
  exact subst_iteratedLubinTate_of_semilinear_intertwine _ (constantCoeff_map_eq_zero _ hf0)
    (constantCoeff_map_eq_zero _ hg0) hθ0 (map_baseIntHom_fixed K σ f)
    (map_baseIntHom_fixed K σ g) hint n

/-! ## 2. `𝒪_{ℂ_K}` の中での評価 —— M3 の器の上での基本性質 -/

/-- 定数項 `0` の冪級数を `𝒪_{ℂ_K}` の位相的冪零元で評価してもノルムは増えない
(`X` をくくり出して `‖z·h(z)‖ = ‖z‖·‖h(z)‖ ≤ ‖z‖`)。 -/
theorem norm_aeval_closureCompletionInt_le (K : PAdicLocalField p)
    {z : ↥(closureCompletionInt K)} (hz : PowerSeries.HasEval z)
    (θ : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0) :
    ‖PowerSeries.aeval hz θ‖ ≤ ‖z‖ := by
  obtain ⟨h, hh⟩ := PowerSeries.X_dvd_iff.mpr hθ0
  rw [hh, map_mul, aeval_X_eq_self]
  rw [norm_coe_closureCompletionInt, norm_coe_closureCompletionInt, Subring.coe_mul, norm_mul]
  calc ‖(z : closureCompletion K)‖
        * ‖((PowerSeries.aeval hz h : ↥(closureCompletionInt K)) : closureCompletion K)‖
      ≤ ‖(z : closureCompletion K)‖ * 1 :=
        mul_le_mul_of_nonneg_left (norm_le_one_closureCompletionInt K _) (norm_nonneg _)
    _ = _ := mul_one _

/-- `0 ∈ 𝒪_{ℂ_K}` は位相的冪零。 -/
theorem hasEval_zero_closureCompletionInt (K : PAdicLocalField p) :
    PowerSeries.HasEval (0 : ↥(closureCompletionInt K)) :=
  tendsto_pow_atTop_nhds_zero_of_norm_lt_one (by rw [norm_zero]; exact zero_lt_one)

/-- 評価点を等しいものへ取り替えてよい(`PowerSeries.HasEval` は `Prop` なので
証明の取り方に依らない)。 -/
theorem aeval_congr_point (K : PAdicLocalField p) {z w : ↥(closureCompletionInt K)}
    (hz : PowerSeries.HasEval z) (hw : PowerSeries.HasEval w) (h : z = w)
    (F : PowerSeries ↥(unramifiedCompletionInt K)) :
    PowerSeries.aeval hz F = PowerSeries.aeval hw F := by
  subst h; rfl

/-- 定数項 `0` の冪級数の `0` での値は `0`。 -/
theorem aeval_zero_eq_zero (K : PAdicLocalField p)
    (h0 : PowerSeries.HasEval (0 : ↥(closureCompletionInt K)))
    (F : PowerSeries ↥(unramifiedCompletionInt K)) (hF : PowerSeries.constantCoeff F = 0) :
    PowerSeries.aeval h0 F = 0 := by
  obtain ⟨h, hh⟩ := PowerSeries.X_dvd_iff.mpr hF
  rw [hh, map_mul, aeval_X_eq_self, zero_mul]

/-- **`‖θ(λ)‖ ≤ ‖λ‖ < 1`**。 -/
theorem norm_evalAt_lt_one (K : PAdicLocalField p) {lam : K.closure} (hlam : ‖lam‖ < 1)
    (θ : PowerSeries ↥(unramifiedCompletionInt K)) (hθ0 : PowerSeries.constantCoeff θ = 0) :
    ‖evalAt K θ hlam‖ < 1 := by
  refine lt_of_le_of_lt (norm_aeval_closureCompletionInt_le K
    (hasEval_coe_of_norm_lt_one K hlam) θ hθ0) ?_
  rw [norm_coe_closureCompletionInt]
  show ‖(lam : closureCompletion K)‖ < 1
  rwa [norm_coe_closureCompletion]

/-- `θ(λ)` はまた位相的冪零——これで `θ(λ)` を新しい評価点として使える。 -/
theorem hasEval_evalAt (K : PAdicLocalField p) {lam : K.closure} (hlam : ‖lam‖ < 1)
    (θ : PowerSeries ↥(unramifiedCompletionInt K)) (hθ0 : PowerSeries.constantCoeff θ = 0) :
    PowerSeries.HasEval (evalAt K θ hlam) :=
  tendsto_pow_atTop_nhds_zero_of_norm_lt_one (norm_evalAt_lt_one K hlam θ hθ0)

/-! ## 3. `𝒪_K` 係数の多項式を `𝒪_{ℂ_K}` で評価する -/

/-- `𝒪_K → 𝒪_{ℂ_K}` は `𝒪_K → K^al → ℂ_K` と一致する。 -/
theorem coe_algebraMap_base_closureCompletionInt (K : PAdicLocalField p) (a : 𝒪[K.carrier]) :
    ((algebraMap (𝒪[K.carrier]) ↥(closureCompletionInt K) a : ↥(closureCompletionInt K))
        : closureCompletion K)
      = closureCompletionCoe K (algebraMap (𝒪[K.carrier]) K.closure a) := by
  rw [algebraMap_base_closureCompletionInt_coe, UniformSpace.Completion.algebraMap_def]
  rfl

/-- `𝒪_K` 係数の多項式を `𝒪_{K̂^{ur}}` 係数の冪級数と見て評価しても、
`𝒪_K`-代数としての多項式評価と同じ(スカラー塔 `𝒪_K → 𝒪_{K̂^{ur}} → 𝒪_{ℂ_K}`)。 -/
theorem aeval_map_baseIntHom_coe (K : PAdicLocalField p) (P : Polynomial (𝒪[K.carrier]))
    {z : ↥(closureCompletionInt K)} (hz : PowerSeries.HasEval z) :
    PowerSeries.aeval hz (PowerSeries.map (baseIntHom K) (P : PowerSeries (𝒪[K.carrier])))
      = Polynomial.aeval z P := by
  rw [← Polynomial.polynomial_map_coe, PowerSeries.aeval_coe]
  show Polynomial.aeval z (Polynomial.map
    (algebraMap (𝒪[K.carrier]) ↥(unramifiedCompletionInt K)) P) = _
  rw [Polynomial.aeval_map_algebraMap]

/-- 評価点が `λ ∈ K^al` の像であるとき、`𝒪_K` 係数多項式の値は `P(λ)` の像。 -/
theorem coe_aeval_closureCompletionInt (K : PAdicLocalField p)
    (P : Polynomial (𝒪[K.carrier])) {z : ↥(closureCompletionInt K)} {lam : K.closure}
    (hz : (z : closureCompletion K) = closureCompletionCoe K lam) :
    ((Polynomial.aeval z P : ↥(closureCompletionInt K)) : closureCompletion K)
      = closureCompletionCoe K ((P.map (algebraMap (𝒪[K.carrier]) K.closure)).eval lam) := by
  rw [Polynomial.eval_map, Polynomial.aeval_def]
  have h1 := Polynomial.hom_eval₂ P (algebraMap (𝒪[K.carrier]) ↥(closureCompletionInt K))
    (SubringClass.subtype (closureCompletionInt K)) z
  have h2 := Polynomial.hom_eval₂ P (algebraMap (𝒪[K.carrier]) K.closure)
    (closureCompletionCoe K) lam
  rw [show (SubringClass.subtype (closureCompletionInt K)).comp
      (algebraMap (𝒪[K.carrier]) ↥(closureCompletionInt K))
      = (closureCompletionCoe K).comp (algebraMap (𝒪[K.carrier]) K.closure) from
    RingHom.ext fun a => coe_algebraMap_base_closureCompletionInt K a,
    show (SubringClass.subtype (closureCompletionInt K)) z = closureCompletionCoe K lam
      from hz] at h1
  exact h1.trans h2.symm

/-- 一般の評価点についての同じ書き換え(`ℂ_K` の中の多項式の値として読む)。 -/
theorem coe_aeval_closureCompletionInt_eq_eval (K : PAdicLocalField p)
    (P : Polynomial (𝒪[K.carrier])) (z : ↥(closureCompletionInt K)) :
    ((Polynomial.aeval z P : ↥(closureCompletionInt K)) : closureCompletion K)
      = Polynomial.eval (z : closureCompletion K)
          (Polynomial.map (closureCompletionCoe K)
            (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure) P)) := by
  rw [Polynomial.map_map, Polynomial.eval_map, Polynomial.aeval_def]
  have h1 := Polynomial.hom_eval₂ P (algebraMap (𝒪[K.carrier]) ↥(closureCompletionInt K))
    (SubringClass.subtype (closureCompletionInt K)) z
  rw [show (SubringClass.subtype (closureCompletionInt K)).comp
      (algebraMap (𝒪[K.carrier]) ↥(closureCompletionInt K))
      = (closureCompletionCoe K).comp (algebraMap (𝒪[K.carrier]) K.closure) from
    RingHom.ext fun a => coe_algebraMap_base_closureCompletionInt K a] at h1
  exact h1

/-! ## 4. `[π^n]_f(z) = 0` と `D_n(z) = 0` の往復、そして根の同定 -/

/-- **`[π^n]_f(z) = 0 ↔ D_n(z) = 0`**(`𝒪_{ℂ_K}` 版)。
Weierstrass 分解 `f^{(n)} = D_n · U_n` の単数部分 `U_n` を約すだけ
(`AdjoinIntegers.lean::eq_zero_of_pi_pow_action_eq_zero` の `𝒪_{ℂ_K}` 版で、
逆向きも同時に取れる形)。 -/
theorem aeval_map_baseIntHom_iteratedLubinTate_eq_zero_iff (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) {z : ↥(closureCompletionInt K)} (hz : PowerSeries.HasEval z) :
    PowerSeries.aeval hz (PowerSeries.map (baseIntHom K) (iteratedLubinTate f n)) = 0
      ↔ Polynomial.aeval z (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n) = 0 := by
  rw [iteratedLubinTate_eq_distinguished_mul_unit hq hπmax hπne0 f hf0 hf1 hf n, map_mul, map_mul,
    IsUnit.mul_left_eq_zero, aeval_map_baseIntHom_coe]
  exact ((isUnit_iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf n).map
    (PowerSeries.map (baseIntHom K))).map (PowerSeries.aeval hz)

/-- `λ ∈ Λ_{f,n}` ならば `f^{(n)}(λ) = 0`(`𝒪_{ℂ_K}` の中で)。 -/
theorem aeval_map_baseIntHom_iteratedLubinTate_eq_zero (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) {z : ↥(closureCompletionInt K)} (hz : PowerSeries.HasEval z) {lam : K.closure}
    (hzlam : (z : closureCompletion K) = closureCompletionCoe K lam)
    (hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    PowerSeries.aeval hz (PowerSeries.map (baseIntHom K) (iteratedLubinTate f n)) = 0 := by
  have hroot : (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).eval lam = 0 := by
    rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots'] at hlam
    exact hlam.2
  rw [aeval_map_baseIntHom_iteratedLubinTate_eq_zero_iff K hq hπmax hπne0 f hf0 hf1 hf n hz]
  apply Subtype.ext
  rw [coe_aeval_closureCompletionInt K _ hzlam, hroot, map_zero]
  rfl

/-- ★★★★★★★★★★★★**代数閉性を使わない根の同定**——`D_n(y) = 0`
(`y ∈ 𝒪_{ℂ_K}`)ならば `y` は `Λ_n` の元の像。

抽象核 `exists_mem_roots_of_eval_map_eq_zero` に
`Q := D_n^{K^al}`・`g := ι_al` を代入するだけ。`hcard` は
`card_roots_iteratedLubinTateDistinguished_map`(根はちょうど `q^n` 個)が
供給する。 -/
theorem exists_mem_torsionPoints_of_aeval_eq_zero (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) {y : ↥(closureCompletionInt K)}
    (hy0 : Polynomial.aeval y (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n) = 0) :
    ∃ r ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n,
      (y : closureCompletion K) = closureCompletionCoe K r := by
  have hmonicD :=
    (isDistinguishedAt_iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n).monic
  have hmonic : (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).Monic := hmonicD.map _
  have hcard : Multiset.card (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).roots
      = (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).natDegree := by
    rw [card_roots_iteratedLubinTateDistinguished_map K hq hπmax hπne0 f hf0 hf1 hf n,
      hmonicD.natDegree_map,
      natDegree_iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n]
  have heval : (Polynomial.map (closureCompletionCoe K)
      (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
        (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n))).eval
      (y : closureCompletion K) = 0 := by
    rw [← coe_aeval_closureCompletionInt_eq_eval, hy0]
    rfl
  obtain ⟨r, hr, hry⟩ :=
    exists_mem_roots_of_eval_map_eq_zero (closureCompletionCoe K) hmonic hcard heval
  exact ⟨r, Multiset.mem_toFinset.mpr hr, hry⟩

/-! ## 5. #7b —— `f_n(λ) = 0 ⟹ g_n(θλ) = 0 ⟹ θλ ∈ ι(Λ_{g,n})` -/

/-- ★★★★★★★★★★★★★★★★★★**#7b**——`θ` は `Λ_{f,n}` を `Λ_{g,n}` へ送る。

原典 51 ページの
「for any α ∈ K^al, f_n(α) = 0 ⟹ g_n(θ(α)) = 0」そのもの。段取り:

1. #7a で `(σ^nθ) ∘ f^{(n)} = g^{(n)} ∘ θ`。
2. 両辺を `λ` で評価する(`aeval_subst_eq_aeval_aeval` = 連鎖律)。
   左辺は `f^{(n)}(λ) = 0` なので `(σ^nθ)(0) = 0`
   (`aeval_zero_eq_zero`、ここで `θ` の定数項が `0` であることを使う)。
3. よって `g^{(n)}(θ(λ)) = 0`。単数部分を約して `D^g_n(θ(λ)) = 0`。
4. `exists_mem_torsionPoints_of_aeval_eq_zero`(整域性だけ)で
   `θ(λ) = ι r`、`r ∈ Λ_{g,n}`。

★`f` と `g` は別々の素元 `π`・`ϖ` に属してよい(`hπmax`・`hϖmax` は別引数)。
`hint` が両者を繋いでいる。 -/
theorem exists_mem_torsionPoints_evalAt (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (σ : unramGal K) (θ : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hint : PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
        (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g))
    (n : ℕ) (lam : K.closure)
    (hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    ∃ r ∈ iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n,
      ((evalAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n θ lam hlam : ↥(closureCompletionInt K))
        : closureCompletion K) = closureCompletionCoe K r := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hg0c : PowerSeries.constantCoeff g = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hg0
  have hnorm := norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n
    lam hlam
  have hz : PowerSeries.HasEval
      (⟨(lam : closureCompletion K), (mem_closureCompletionInt K _).mpr
        (by rw [norm_coe_closureCompletion]; exact hnorm.le)⟩ : ↥(closureCompletionInt K)) :=
    hasEval_coe_of_norm_lt_one K hnorm
  have hy : PowerSeries.HasEval (evalAt K θ hnorm) := hasEval_evalAt K hnorm θ hθ0
  have hw0 : PowerSeries.aeval hz (PowerSeries.map (baseIntHom K) (iteratedLubinTate f n)) = 0 :=
    aeval_map_baseIntHom_iteratedLubinTate_eq_zero K hq hπmax hπne0 f hf0 hf1 hf n hz rfl hlam
  have hkey := congrArg (PowerSeries.aeval hz)
    (subst_iteratedLubinTate_of_dworkTheta K σ hf0c hg0c hθ0 hint n)
  rw [aeval_subst_eq_aeval_aeval
      (PowerSeries.HasSubst.of_constantCoeff_zero'
        (constantCoeff_map_eq_zero _
          (constantCoeff_iteratedLubinTate_of_constantCoeff_zero hf0c n)))
      (constantCoeff_map_eq_zero _
        (constantCoeff_iteratedLubinTate_of_constantCoeff_zero hf0c n))
      hz (hasEval_zero_closureCompletionInt K) hw0,
    aeval_zero_eq_zero K _ _ (constantCoeff_iterate_map_eq_zero _ hθ0 n),
    aeval_subst_eq_aeval_aeval (PowerSeries.HasSubst.of_constantCoeff_zero' hθ0) hθ0 hz hy rfl]
    at hkey
  have hD : Polynomial.aeval (evalAt K θ hnorm)
      (iteratedLubinTateDistinguished hq hϖmax hϖne0 g hg0 hg1 hg n) = 0 :=
    (aeval_map_baseIntHom_iteratedLubinTate_eq_zero_iff K hq hϖmax hϖne0 g hg0 hg1 hg n hy).mp
      hkey.symm
  exact exists_mem_torsionPoints_of_aeval_eq_zero K hq hϖmax hϖne0 g hg0 hg1 hg n hD

/-- ★原典の「and, similarly, `g_n(α) = 0 ⟹ f_n(θ^{-1}(α)) = 0`」——
合成逆 `θ'` は `Λ_{g,n}` を `Λ_{f,n}` へ送る。

抽象核 `subst_intertwine_of_comp_inverse` で絡みを裏返し、#7b に
`(f, g, θ) := (g, f, θ')` を代入するだけ。 -/
theorem exists_mem_torsionPoints_evalAt_inverse (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (σ : unramGal K) (θ θ' : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0) (hθ'0 : PowerSeries.constantCoeff θ' = 0)
    (hθθ' : PowerSeries.subst θ' θ = PowerSeries.X)
    (hθ'θ : PowerSeries.subst θ θ' = PowerSeries.X)
    (hint : PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
        (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g))
    (n : ℕ) (mu : K.closure)
    (hmu : mu ∈ iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n) :
    ∃ r ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n,
      ((evalAtTorsionPoint K hq hϖmax hϖne0 g hg0 hg1 hg n θ' mu hmu
        : ↥(closureCompletionInt K)) : closureCompletion K) = closureCompletionCoe K r := by
  refine exists_mem_torsionPoints_evalAt K hq hϖmax hϖne0 g hg0 hg1 hg hπmax hπne0 f hf0 hf1 hf
    σ θ' hθ'0 ?_ n mu hmu
  exact subst_intertwine_of_comp_inverse _
    (constantCoeff_map_eq_zero _
      (by rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0))
    hθ0 hθ'0 hθθ' hθ'θ hint

/-! ## 6. #7c —— 全単射 -/

/-- **`θ'(θ(λ)) = λ`** —— 合成逆の評価。連鎖律 `aeval_subst_eq_aeval_aeval` と
`subst θ θ' = X` から直ちに従う。★これが `θ` の単射性の全て。 -/
theorem coe_aeval_evalAt_comp_inverse (K : PAdicLocalField p) {lam : K.closure}
    (hnorm : ‖lam‖ < 1) (θ θ' : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hθ'θ : PowerSeries.subst θ θ' = PowerSeries.X) :
    ((PowerSeries.aeval (hasEval_evalAt K hnorm θ hθ0) θ' : ↥(closureCompletionInt K))
      : closureCompletion K) = closureCompletionCoe K lam := by
  have h := aeval_subst_eq_aeval_aeval (p := θ') (PowerSeries.HasSubst.of_constantCoeff_zero' hθ0)
    hθ0 (hasEval_coe_of_norm_lt_one K hnorm) (hasEval_evalAt K hnorm θ hθ0) rfl
  rw [hθ'θ, aeval_X_eq_self] at h
  rw [← h]
  rfl

/-- ★★★★★★★★★★★★★★★★★★★★**#7c**——`θ` は
`Λ_{f,n}` から `Λ_{g,n}` への**全単射**を誘導する。

原典 51 ページ「Therefore θ defines a bijection Λ_{f,n} → Λ_{g,n}」。

* **写る**(`MapsTo`): #7b。
* **単射**(`InjOn`): `θ'(θ(λ)) = λ`(`coe_aeval_evalAt_comp_inverse`)と
  `ι_al` の単射性。
* **全射**(`SurjOn`): `|Λ_{f,n}| = |Λ_{g,n}| = q^n`
  (`card_iteratedLubinTateTorsionPoints`)と単射性から
  `Finset.eq_of_subset_of_card_le`。

★★逸脱(記録): 原典は全射を「逆向きの包含」から出すが、ここでは
濃度から出した。逆向きの包含も `exists_mem_torsionPoints_evalAt_inverse`
として別に証明してある。 -/
theorem exists_bijOn_torsionPoints_of_dworkTheta (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (σ : unramGal K) (θ θ' : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hθ'θ : PowerSeries.subst θ θ' = PowerSeries.X)
    (hint : PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
        (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g))
    (n : ℕ) :
    ∃ e : K.closure → K.closure,
      (∀ (lam : K.closure)
          (hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n),
        ((evalAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n θ lam hlam
            : ↥(closureCompletionInt K)) : closureCompletion K)
          = closureCompletionCoe K (e lam)) ∧
      Set.BijOn e ↑(iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
        ↑(iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n) := by
  classical
  have H : ∀ lam : K.closure, ∃ r : K.closure,
      lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n →
        (r ∈ iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n ∧
          ∀ hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n,
            ((evalAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n θ lam hlam
              : ↥(closureCompletionInt K)) : closureCompletion K) = closureCompletionCoe K r) := by
    intro lam
    by_cases hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n
    · obtain ⟨r, hr, hrv⟩ := exists_mem_torsionPoints_evalAt K hq hπmax hπne0 f hf0 hf1 hf
        hϖmax hϖne0 g hg0 hg1 hg σ θ hθ0 hint n lam hlam
      exact ⟨r, fun _ => ⟨hr, fun _ => hrv⟩⟩
    · exact ⟨0, fun h => absurd h hlam⟩
  choose e he using H
  have hmaps : Set.MapsTo e ↑(iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
      ↑(iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n) :=
    fun lam hlam => (he lam (Finset.mem_coe.mp hlam)).1
  have hinj : Set.InjOn e ↑(iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) := by
    intro a ha b hb hab
    have ha' := Finset.mem_coe.mp ha
    have hb' := Finset.mem_coe.mp hb
    have hna := norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n
      a ha'
    have hnb := norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n
      b hb'
    have hva : ((evalAt K θ hna : ↥(closureCompletionInt K)) : closureCompletion K)
        = closureCompletionCoe K (e a) := (he a ha').2 ha'
    have hvb : ((evalAt K θ hnb : ↥(closureCompletionInt K)) : closureCompletion K)
        = closureCompletionCoe K (e b) := (he b hb').2 hb'
    have hy : evalAt K θ hna = evalAt K θ hnb := Subtype.ext (by rw [hva, hvb, hab])
    apply closureCompletionCoe_injective K
    rw [← coe_aeval_evalAt_comp_inverse K hna θ θ' hθ0 hθ'θ,
      ← coe_aeval_evalAt_comp_inverse K hnb θ θ' hθ0 hθ'θ,
      aeval_congr_point K (hasEval_evalAt K hna θ hθ0) (hasEval_evalAt K hnb θ hθ0) hy θ']
  have himg : (iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n).image e
      = iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n := by
    refine Finset.eq_of_subset_of_card_le ?_ ?_
    · intro x hx
      obtain ⟨a, ha, rfl⟩ := Finset.mem_image.mp hx
      exact (he a ha).1
    · rw [Finset.card_image_of_injOn hinj,
        card_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n,
        card_iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n]
  refine ⟨e, fun lam hlam => (he lam hlam).2 hlam, hmaps, hinj, ?_⟩
  show (↑(iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n) : Set K.closure)
    ⊆ e '' ↑(iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
  rw [← himg, Finset.coe_image]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**#7 の完成形**——`θ` を仮定に置かず、
Λ6(Dwork)から取ってくる形。

`f ∈ F_π` と `g ∈ F_ϖ`(`ϖ = uπ`、`u` は単数)について、
`Λ_{f,n}` と `Λ_{g,n}` は**全単射**である。

★退化の自己検査: `u` が単数であることを落とすと `span {uπ}` が
`span {π}` と一致せず、そもそも `Λ_{g,n}` の定義に必要な
`maximalIdeal = span {uπ}` が言えない。すなわち「同じ素元の塔」という
条件は statement に埋め込まれている。 -/
theorem exists_bijOn_iteratedLubinTateTorsionPoints (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {u : 𝒪[K.carrier]} (hu : IsUnit u)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = u * π)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    ∃ e : K.closure → K.closure,
      Set.BijOn e ↑(iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
        ↑(iteratedLubinTateTorsionPoints K hq
            (hπmax.trans (Ideal.span_singleton_mul_left_unit hu π).symm)
            (mul_ne_zero hu.ne_zero hπne0) g hg0 hg1 hg n) := by
  obtain ⟨σ, -, -, -, hstep2⟩ := exists_arithFrobenius_isCoherent_dworkThetaStep2 K hq
  obtain ⟨θ, θ', hθ0, -, -, hθ'θ, -, hint⟩ :=
    hstep2 π hπmax f hf0 hf1 hf u hu g hg0 hg1 hg
  obtain ⟨e, -, hbij⟩ := exists_bijOn_torsionPoints_of_dworkTheta K hq hπmax hπne0 f hf0 hf1 hf
    (hπmax.trans (Ideal.span_singleton_mul_left_unit hu π).symm)
    (mul_ne_zero hu.ne_zero hπne0) g hg0 hg1 hg σ θ θ'
    (by rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hθ0) hθ'θ hint n
  exact ⟨e, hbij⟩

end ABC3.Found.PGC
