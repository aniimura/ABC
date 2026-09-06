import ABC3.Found.PGC.DworkMultiplicative

/-!
# Dwork の補題(核)—— `σ` の固定環はちょうど `𝒪_K` である(`sorry` 無し)

経路 Λ の節点 **Λ6b(M2)**。`DworkAdditive.lean` / `DworkMultiplicative.lean`
(第 1059)が出したのは Dwork の補題の**全射性**だけであった
(`σb - b = c` と `σξ/ξ = u` の可解性)。原典 LEMMA 3.11 のもう半分、すなわち
**核の記述**

```
∃ σ ∈ Gal(K^ur/K), (σ は剰余体で z ↦ z^q) ∧ {b ∈ 𝒪_{K̂^ur} | σb = b} = 𝒪_K の像
```

が本ファイルである。Milne CFT の **Prop 3.10 Step 2**(「`σh = h` ゆえ
`h ∈ 𝒪_K[[T]]`」)が使うのはこちらの向きなので、Λ6 の Step 2 の入口にあたる。

★**量化の順が主張の中身**である。`∃ σ, (剰余体で q 乗) ∧ (固定環の記述)` を
主張しており、`σ` は `∃` の内側に閉じ込めてある。結論に自由なパラメータは出ない
(一意化元 `π` も証明の内側にある)。

## 到達点

| | 主張 |
|---|---|
| 0 | `baseIntHom`:`𝒪_K → 𝒪_{K̂^{ur}}`(単射・局所準同型) |
| 1 | `exists_residue_baseIntHom_eq`:**剰余体の段**——`t^q = t` なら `t` は `𝓀_K` の像 |
| 2 | `exists_sub_baseIntHom_eq_uniformizer_mul`:1 段の補正(`b - a₀ = π b₁`、`b₁` も固定) |
| 3 | `exists_sub_baseIntHom_eq_uniformizer_pow_mul`:`N` 段の近似(`N` の帰納法 1 本) |
| 4 | `mem_span_pow_of_baseIntHom_mem`:**降下**——`𝒪_{K̂^{ur}}` で `(π^N)` なら `𝒪_K` でも |
| 5 | `exists_baseIntHom_eq_of_unramGalCompletionInt_eq`:**固定 ⇒ `𝒪_K` の像** |
| 6 | `unramGalCompletionInt_eq_self_iff` / `setOf_unramGalCompletionInt_eq_self`:環版の同値 |
| 7 | `unramGalCompletionUnits_eq_self_iff` / `setOf_unramGalCompletionUnits_eq_self`:**単元版** |
| 8 | `eqLocus_unramGalCompletionInt`:`Subring` の等式(`RingHom.eqLocus`) |
| 9 | **`exists_arithFrobenius_dworkFixedRing`** / `exists_arithFrobenius_dworkFixedUnits` |
| 10 | `exists_arithFrobenius_dworkKernelAndSurjective`:核と全射性を**同じ 1 つの `σ`** で |

## ★見立てとの差(2026-09-06 の実測)

段取り係の見立ては 2 つあった。**両方とも当たっている**。

* **「剰余体の段は `mem_of_pow_eq_self_of_card_eq` で済む」——当たり。**
  ただし在庫の補題は「`𝓀_{K̂^{ur}}` の Frobenius 固定部分 = `𝓀_K`」を
  **直接には言っていない**(中身は「`X^m - X` の根はたかだか `m` 個」という
  純粋に多項式の補題である)。当てるには `T := 𝓀_K の像`(ちょうど `q` 元)を
  自分で作る必要があった。それでも `exists_residue_baseIntHom_eq` は **20 行**で、
  見積が下がったという判断は正しい。
* **「`𝒪_K` の像が閉じている(`IsAdicComplete`)の移送が最初の関門かもしれない」
  ——外れ(思ったより安い)。** `𝒪_K` の `𝔪` 進完備性は加法版と**同じ 3 行**で出た:

  ```lean
  haveI : IsDiscreteValuationRing ((Valued.v : Valuation K.carrier NNReal).valuationSubring) :=
    isDiscreteValuationRing_carrierIntegers K
  exact isAdicComplete_valuationSubring
  ```

  `𝒪[K.carrier]` は `Valued.integer K.carrier`(`Subring`)、
  `isAdicComplete_valuationSubring` は `Valued.v.valuationSubring`(`ValuationSubring`)
  についてだが、**既定透明度で defeq** なので `exact` がそのまま通る。
  ★本当に効いたのは完備性の移送ではなく**降下**(段 4)の方で、そこは
  `IsAdicComplete` ではなく**ノルム 1 本**(`‖c‖ ≤ ‖π‖^N` から `c/π^N ∈ 𝒪_K`)で済んだ。

全体で `lean_check` の往復は **10 回**、うち失敗は **2 回**(どちらも配管:
`rw ... at *` の誤爆と、部分環の `↑(a-b)` を `sub_eq_zero` に食わせた形)。

## 筋(6 段)

1. `baseIntHom K : 𝒪_K →+* 𝒪_{K̂^{ur}}`(`algebraMap` の制限)。単射・`IsLocalHom`。
2. **剰余体**:`σb = b` なら `residue b` は `q` 乗で不変。`𝓀_K` の像は
   ちょうど `q` 元で全部 `q` 乗不変なので、`mem_of_pow_eq_self_of_card_eq`
   (`X^q - X` の根はたかだか `q` 個)で `residue b` は `𝓀_K` の像に入る。
3. **1 段**:代表 `a₀ ∈ 𝒪_K` を取ると `b - a₀ ∈ 𝔪 = (π)`、すなわち `b - a₀ = π b₁`。
   `σπ = π` と `σa₀ = a₀` から `π(σb₁ - b₁) = 0`、`𝒪_{K̂^{ur}}` は整域なので
   **`b₁` も固定される**。★ここが再帰を回す鍵。
4. **`N` 段**:`b` を全称量化したまま `N` の帰納法。`a := a₀ + π a₁` と組み替えるだけ
   (`linear_combination` 1 行)。`choose` も `Function.iterate` も要らない。
5. **降下**:`baseIntHom c ∈ (π^N)` なら `‖c‖ ≤ ‖π‖^N` なので `c/π^N ∈ 𝒪_K`、
   すなわち `c ∈ (π^N)`。★これで近似列 `a_N` が `𝒪_K` の中で Cauchy になる。
6. **極限**:`𝒪_K` の `𝔪` 進完備性で `a := lim a_N` を取り、
   `b - a ∈ (π^N)`(∀N)を `𝒪_{K̂^{ur}}` の `IsHausdorff` で `0` にする。

## 材料(すべて本プロジェクトの在庫)

| 要るもの | 出どころ |
|---|---|
| `X^m - X` の根はたかだか `m` 個 | `UnramifiedResidueField.mem_of_pow_eq_self_of_card_eq` |
| `σ` の剰余体での作用 | `DworkAdditive.residue_unramGalCompletionInt` |
| `σπ = π` | `DworkAdditive.unramGalCompletionInt_uniformizerCompletionInt` |
| `𝔪 = (π)` / `𝔪^N = (π^N)` | `UnramifiedCompletionDVR` / `DworkAdditive` |
| `𝒪_{K̂^{ur}}` の `𝔪` 進完備性 | `DworkAdditive.isAdicComplete_unramifiedCompletionInt` |
| `𝒪_K` の DVR 性 | `UnramifiedExtension.isDiscreteValuationRing_carrierIntegers` |
| `𝔪` 進完備性(一般形) | `GaloisRep.isAdicComplete_valuationSubring` |
| 算術 Frobenius | `ArithFrobeniusTopGen.arithFrobenius` / `arithFrobenius_residue` |
| 全射性(相方) | `DworkAdditive` / `DworkMultiplicative` |

★**`frobeniusCompletionInt` は使っていない**。あれは `coherentFrobenius`
(位相的生成元一般)から作られており、剰余体上の作用は `z ↦ z^{q^k}` でよい。
固定環の主張は**算術 Frobenius でなければ変わる**
(`UnramifiedResidueField.lean` の訂正欄)。

## ★設計上の注意(守ったこと)

* **既存の `Found/PGC/*.lean` を書き換えていない**(import のみ)。
* **`Skeleton` / `Interface` を触っていない**。`inertia` を経由せず、
  `Interface` の `SubgroupCorrespondence` / `ResidueCardinality` は
  主張にも証明にも現れない(Corollary 1.3 との循環を避けるため)。
* **結論に自由なパラメータを出していない**——主定理の型は `K` にしか依存せず、
  `σ` も一意化元 `π` も `∃` / 証明の内側にある。
* **`sorry` 本体の `def` を作っていない**(D21)。

## 退化の自己検査

* **`σ` を任意の `unramGal K` に替えると偽**。`σ = 1` なら固定環は `𝒪_{K̂^{ur}}`
  全体で、`𝒪_K` の像より真に大きい(`𝒪_{K̂^{ur}}` は `𝒪_K` 上無限次)。
  ★だから `σ` は `∃` の内側に無ければならない。
* **`σ` を「位相的生成元」に弱めた場合に何が起きるかは未確認**。
  この証明は回らない(段 2 の剰余方程式が `t^{q^k} = t`(`k ∈ Ẑ`)になり
  多項式方程式でなくなるので `mem_of_pow_eq_self_of_card_eq` を当てられない)が、
  ★**主張自体が偽になるかは確かめていない**。「偽」とは書かない。
* `hσ`(剰余体で `q` 乗)を落とすと段 2 が崩れる。`⊇`(`𝒪_K` の像は固定される)は
  `hσ` 無しでも真だが、`⊆` は上のとおり `σ = 1` で反例になる。
* `𝒪_{K̂^{ur}}` を `K̂^{ur}` に替えると固定体は `K` そのものになる(真だと思われるが
  下流は `𝒪` 版しか使わないので示していない)。
* `K^{ur}` を `K.closure` に替えると段 6 が**偽**——`𝒪_{ℂ_p}` は DVR でなく
  `𝔪` 進完備でもない。効いているのは不分岐性である。
* 単元版の左辺を「固定される `0` でない元」に緩めると**偽**——`π` は固定される
  `0` でない元だが `𝒪_K^×` の像には入らない(環版の主張なら真)。

## 逸脱(記録)

* **原典(Milne, CFT の LEMMA 3.11)の証明路を取っていない**。原典は
  `𝒪_{K̂^{ur}}` を逆極限として作り、各段 `𝒪/𝔪^n` で完全列を回して
  A.7 / A.8(逆極限の完全性)に渡す。本ファイルはその路を**使わない**
  ——Λ5 の逸脱(「ノルムの完備化と `𝔪` 進完備化の一致は未証明」)がまだ
  閉じていないためである。代わりに `IsAdicComplete` を位相の側から取り、
  逐次近似はイデアルの言葉だけで回した。★得られる主張は同じで、下流への影響は無い。
* Λ5 の逸脱そのもの(ノルム完備化 ≅ `𝔪` 進完備化)は**依然として閉じていない**。
  本ファイルはそれを**必要としない**形に組み替えただけである。
* 段 5(降下)だけは**ノルムを使っている**(加法版・乗法版は一度もノルムを
  使わなかった)。`𝒪_K → 𝒪_{K̂^{ur}}` が忠実平坦であることを代数的に出すより、
  等長性 `‖baseIntHom a‖ = ‖a‖` から 3 行で出す方が安いと判断した。
  ★配管上の選択であり、数学は同じである。
* **`𝒪_K` の像を `Subring` としてではなく `RingHom.range` として書いている**。
  `Set` 版(`setOf_…`)と `Subring` 版(`eqLocus_…`)の両方を出してあるので、
  下流はどちらでも受け取れる。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found.GaloisRep
open scoped NNReal Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 0. `𝒪_K → 𝒪_{K̂^{ur}}`

固定環の「あるべき姿」を名前で 1 度だけ固定する。`uniformizerCompletionInt` は
同じ写像の値だけを見た版なので、両者が `rfl` で繋がることも記録しておく
(`baseIntHom_eq_uniformizerCompletionInt`)。 -/

/-- **`𝒪_K → 𝒪_{K̂^{ur}}`**——`algebraMap K.carrier (K̂^{ur})` の整数環への制限。

★固定環の記述はこの写像の像として述べる。 -/
noncomputable def baseIntHom (K : PAdicLocalField p) :
    𝒪[K.carrier] →+* ↥(unramifiedCompletionInt K) :=
  RingHom.codRestrict
    ((algebraMap K.carrier (unramifiedCompletion K)).comp
      (SubringClass.subtype 𝒪[K.carrier]))
    (unramifiedCompletionInt K)
    (fun x => by
      rw [mem_unramifiedCompletionInt]
      show ‖algebraMap K.carrier (unramifiedCompletion K) (x : K.carrier)‖ ≤ 1
      rw [norm_algebraMap_unramifiedCompletion]
      have h := x.2; rw [Valued.integer.mem_iff] at h; exact h)

@[simp] theorem baseIntHom_coe (K : PAdicLocalField p) (a : 𝒪[K.carrier]) :
    ((baseIntHom K a : ↥(unramifiedCompletionInt K)) : unramifiedCompletion K)
      = algebraMap K.carrier (unramifiedCompletion K) (a : K.carrier) := rfl

/-- `uniformizerCompletionInt`(`UnramifiedCompletionDVR`)は `baseIntHom` の値である。 -/
theorem baseIntHom_eq_uniformizerCompletionInt (K : PAdicLocalField p) (a : 𝒪[K.carrier]) :
    baseIntHom K a = uniformizerCompletionInt K a := rfl

/-- **`baseIntHom` は等長**——`K̂^{ur}` のノルムは `K` のノルムを延長するから。 -/
theorem norm_baseIntHom (K : PAdicLocalField p) (a : 𝒪[K.carrier]) :
    ‖((baseIntHom K a : ↥(unramifiedCompletionInt K)) : unramifiedCompletion K)‖
      = ‖(a : K.carrier)‖ := by
  rw [baseIntHom_coe, norm_algebraMap_unramifiedCompletion]

theorem baseIntHom_injective (K : PAdicLocalField p) : Function.Injective (baseIntHom K) := by
  intro a b hab
  have hn : ‖((baseIntHom K (a - b) : ↥(unramifiedCompletionInt K))
      : unramifiedCompletion K)‖ = 0 := by
    rw [map_sub, hab, sub_self]; simp
  rw [norm_baseIntHom, norm_eq_zero] at hn
  exact sub_eq_zero.mp (Subtype.ext hn)

/-- **`baseIntHom` は局所準同型**——単元版(`§6`)で `a` が `𝒪_K` の単元であることを
言うのに要る。 -/
instance isLocalHom_baseIntHom (K : PAdicLocalField p) : IsLocalHom (baseIntHom K) := by
  constructor
  intro a ha
  by_contra hcon
  rw [Valued.integer.isUnit_iff_norm_eq_one] at hcon
  have hle : ‖(a : K.carrier)‖ ≤ 1 := by
    have h := a.2; rwa [Valued.integer.mem_iff] at h
  refine absurd ha ?_
  rw [not_isUnit_unramifiedCompletionInt, norm_baseIntHom]
  exact lt_of_le_of_ne hle (by simpa using hcon)

/-- **`σ` は `𝒪_K` の像を固定する**(`⊇` の向き)——`σ` は `K`-代数同型だから。 -/
theorem unramGalCompletionInt_baseIntHom (K : PAdicLocalField p) (σ : unramGal K)
    (a : 𝒪[K.carrier]) :
    unramGalCompletionInt K σ (baseIntHom K a) = baseIntHom K a :=
  unramGalCompletionInt_uniformizerCompletionInt K σ a

theorem uniformizerCompletionInt_ne_zero (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπ0 : (π : K.carrier) ≠ 0) : uniformizerCompletionInt K π ≠ 0 := by
  intro h
  apply hπ0
  have hz : ((uniformizerCompletionInt K π : ↥(unramifiedCompletionInt K))
      : unramifiedCompletion K) = 0 := by rw [h]; simp
  rw [uniformizerCompletionInt_coe] at hz
  have hn : ‖algebraMap K.carrier (unramifiedCompletion K) (π : K.carrier)‖ = 0 := by
    rw [hz]; simp
  rw [norm_algebraMap_unramifiedCompletion, norm_eq_zero] at hn
  exact hn

/-! ## 1. 剰余体の段 —— `𝓀_{K̂^{ur}}` の Frobenius 固定部分は `𝓀_K` -/

/-- **★★★★★★★★★★★★(Λ6b)剰余体の段**——`t^q = t` なる
`t ∈ 𝓀_{K̂^{ur}}` は `𝓀_K` の像である。

`T := 𝓀_K の像` はちょうど `q` 元(`𝓀_K → 𝓀_{K̂^{ur}}` は体からの環準同型なので単射)
で、各元は `FiniteField.pow_card` により `q` 乗で不変。したがって
`mem_of_pow_eq_self_of_card_eq`(`X^q - X` の根はたかだか `q` 個)が
`t ∈ T` を与える。

★これが「`K^{ur}/K` の剰余体拡大 `𝔽̄_q/𝔽_q` の Frobenius 固定体は `𝔽_q`」の中身
だが、Galois 理論を一切経由していない(根の数え上げだけ)。

退化の自己検査:`ht`(`q` 乗で不変)を落とすと**偽**。
`𝓀_{K̂^{ur}} = 𝔽̄_q` は `𝓀_K = 𝔽_q` より真に大きい。 -/
theorem exists_residue_baseIntHom_eq (K : PAdicLocalField p)
    {t : IsLocalRing.ResidueField ↥(unramifiedCompletionInt K)}
    (ht : t ^ (Nat.card 𝓀[K.carrier]) = t) :
    ∃ a : 𝒪[K.carrier], IsLocalRing.residue _ (baseIntHom K a) = t := by
  classical
  haveI : Fintype 𝓀[K.carrier] := Fintype.ofFinite _
  set φ := IsLocalRing.ResidueField.map (baseIntHom K) with hφ
  set T : Finset (IsLocalRing.ResidueField ↥(unramifiedCompletionInt K)) :=
    Finset.image φ Finset.univ with hT
  have hTcard : T.card = Nat.card 𝓀[K.carrier] := by
    rw [hT, Finset.card_image_of_injective _ φ.injective, Finset.card_univ,
      Nat.card_eq_fintype_card]
  have hTpow : ∀ s ∈ T, s ^ (Nat.card 𝓀[K.carrier]) = s := by
    intro s hs
    rw [hT, Finset.mem_image] at hs
    obtain ⟨u, -, rfl⟩ := hs
    rw [← map_pow, Nat.card_eq_fintype_card, FiniteField.pow_card]
  have hmem := mem_of_pow_eq_self_of_card_eq (one_lt_residueCard K) hTcard hTpow ht
  rw [hT, Finset.mem_image] at hmem
  obtain ⟨u, -, hu⟩ := hmem
  obtain ⟨a, rfl⟩ := IsLocalRing.residue_surjective u
  exact ⟨a, by rw [← hu, hφ, IsLocalRing.ResidueField.map_residue]⟩

/-! ## 2. 1 段の補正 —— 固定元は固定元に落ちる -/

/-- **★★★★★★★★★★★★★★(Λ6b)1 段の補正**——`σb = b` なら
`b - a₀ = π b₁` なる `a₀ ∈ 𝒪_K` と **`σ` で固定される** `b₁` が取れる。

★「`b₁` も固定される」が再帰を回す鍵である。`σ(π b₁) = π b₁` と `σπ = π` から
`π(σb₁ - b₁) = 0`、`𝒪_{K̂^{ur}}` は整域(体の付値環)なので `σb₁ = b₁`。

退化の自己検査。

* `hσ`(剰余体で `q` 乗)を落とすと段 1 が当たらない。実際 `σ = 1` なら
  すべての `b` が固定されるので `b - a₀ ∈ 𝔪` は一般に偽。
* `hπ0` を落とすと `π = 0` で割り算(`mul_right_cancel₀`)が壊れる。
* `hπmax` を落とすと `π` が自由なパラメータになり `𝔪 = (π)` が偽。 -/
theorem exists_sub_baseIntHom_eq_uniformizer_mul (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    {π : 𝒪[K.carrier]} (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    {b : ↥(unramifiedCompletionInt K)} (hb : unramGalCompletionInt K σ b = b) :
    ∃ a : 𝒪[K.carrier], ∃ b' : ↥(unramifiedCompletionInt K),
      unramGalCompletionInt K σ b' = b' ∧
        b - baseIntHom K a = uniformizerCompletionInt K π * b' := by
  have ht : (IsLocalRing.residue _ b) ^ (Nat.card 𝓀[K.carrier]) = IsLocalRing.residue _ b := by
    rw [← residue_unramGalCompletionInt K hσ b, hb]
  obtain ⟨a, ha⟩ := exists_residue_baseIntHom_eq K ht
  have hmem : b - baseIntHom K a ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K) := by
    rw [← IsLocalRing.residue_eq_zero_iff, map_sub, ha, sub_self]
  rw [maximalIdeal_unramifiedCompletionInt_eq_span K hπ0 hπmax,
    Ideal.mem_span_singleton'] at hmem
  obtain ⟨b', hb'⟩ := hmem
  have hfix : unramGalCompletionInt K σ (b - baseIntHom K a) = b - baseIntHom K a := by
    rw [map_sub, hb, unramGalCompletionInt_baseIntHom]
  rw [← hb', map_mul, unramGalCompletionInt_uniformizerCompletionInt] at hfix
  exact ⟨a, b', mul_right_cancel₀ (uniformizerCompletionInt_ne_zero K hπ0) hfix,
    by rw [← hb', mul_comm]⟩

/-! ## 3. `N` 段の近似 -/

/-- **★★★★★★★★★★★★★★★★(Λ6b)`N` 段の近似**——`σb = b` なら
`b - a = π^N b'` なる `a ∈ 𝒪_K` と固定元 `b'` が取れる。

`b` を全称量化したまま `N` について帰納する。帰納段は `a := a₀ + π a₁` と
組み替えるだけで、`linear_combination` 1 行で閉じる。
★`choose` も `Function.iterate` も要らない(加法版・乗法版との違い)。 -/
theorem exists_sub_baseIntHom_eq_uniformizer_pow_mul (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    {π : 𝒪[K.carrier]} (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) :
    ∀ (N : ℕ) (b : ↥(unramifiedCompletionInt K)), unramGalCompletionInt K σ b = b →
      ∃ a : 𝒪[K.carrier], ∃ b' : ↥(unramifiedCompletionInt K),
        unramGalCompletionInt K σ b' = b' ∧
          b - baseIntHom K a = uniformizerCompletionInt K π ^ N * b' := by
  intro N
  induction N with
  | zero => exact fun b hb => ⟨0, b, hb, by simp⟩
  | succ N ih =>
    intro b hb
    obtain ⟨a₀, b₁, hb₁, h₀⟩ :=
      exists_sub_baseIntHom_eq_uniformizer_mul K hσ hπ0 hπmax hb
    obtain ⟨a₁, b₂, hb₂, h₁⟩ := ih b₁ hb₁
    refine ⟨a₀ + π * a₁, b₂, hb₂, ?_⟩
    have hBπ : baseIntHom K π = uniformizerCompletionInt K π := rfl
    rw [map_add, map_mul, hBπ, pow_succ]
    linear_combination h₀ + uniformizerCompletionInt K π * h₁

/-! ## 4. 降下 —— `𝒪_{K̂^{ur}}` の `(π^N)` は `𝒪_K` の `(π^N)` に落ちる -/

/-- **★★★★★★★★★★(Λ6b)降下**——`baseIntHom c ∈ (π^N)` なら `c ∈ (π^N)`。

`baseIntHom` は等長なので `‖c‖ = ‖w‖·‖π‖^N ≤ ‖π‖^N`、したがって
`c/π^N` のノルムは `1` 以下で `𝒪_K` に入る。

★これが無いと近似列 `a_N` が `𝒪_K` の中で Cauchy にならない
(`𝒪_{K̂^{ur}}` の中で Cauchy なだけでは極限が `𝒪_K` に入らない)。

退化の自己検査:`hπ0` を落とすと `π^N = 0` で `c/π^N` が定義できない。 -/
theorem mem_span_pow_of_baseIntHom_mem (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπ0 : (π : K.carrier) ≠ 0) (N : ℕ) {c : 𝒪[K.carrier]}
    (h : baseIntHom K c ∈ Ideal.span {uniformizerCompletionInt K π ^ N}) :
    c ∈ Ideal.span {π ^ N} := by
  rw [Ideal.mem_span_singleton'] at h ⊢
  obtain ⟨w, hw⟩ := h
  have hcoe : ((w : ↥(unramifiedCompletionInt K)) : unramifiedCompletion K)
      * (algebraMap K.carrier (unramifiedCompletion K) (π : K.carrier)) ^ N
      = algebraMap K.carrier (unramifiedCompletion K) (c : K.carrier) := by
    have h2 := congrArg
      (fun z : ↥(unramifiedCompletionInt K) => (z : unramifiedCompletion K)) hw
    simpa using h2
  have hw1 : ‖((w : ↥(unramifiedCompletionInt K)) : unramifiedCompletion K)‖ ≤ 1 :=
    (mem_unramifiedCompletionInt K _).mp w.2
  have hπpos : 0 < ‖(π : K.carrier)‖ := norm_pos_iff.mpr hπ0
  have hn : ‖(c : K.carrier)‖
      = ‖((w : ↥(unramifiedCompletionInt K)) : unramifiedCompletion K)‖
        * ‖(π : K.carrier)‖ ^ N := by
    rw [← norm_algebraMap_unramifiedCompletion, ← hcoe, norm_mul, norm_pow,
      norm_algebraMap_unramifiedCompletion]
  have hle : ‖(c : K.carrier)‖ ≤ ‖(π : K.carrier)‖ ^ N := by
    rw [hn]
    calc ‖((w : ↥(unramifiedCompletionInt K)) : unramifiedCompletion K)‖
          * ‖(π : K.carrier)‖ ^ N
        ≤ 1 * ‖(π : K.carrier)‖ ^ N :=
          mul_le_mul_of_nonneg_right hw1 (pow_nonneg (norm_nonneg _) N)
      _ = ‖(π : K.carrier)‖ ^ N := one_mul _
  have hd : ‖(c : K.carrier) / (π : K.carrier) ^ N‖ ≤ 1 := by
    rw [norm_div, norm_pow, div_le_one (pow_pos hπpos N)]
    exact hle
  refine ⟨⟨_, Valued.integer.mem_iff.mpr hd⟩, Subtype.ext ?_⟩
  show (c : K.carrier) / (π : K.carrier) ^ N * ((π ^ N : 𝒪[K.carrier]) : K.carrier)
      = (c : K.carrier)
  push_cast
  exact div_mul_cancel₀ _ (pow_ne_zero N hπ0)

/-- `𝒪_K` の `(π^N)` は `𝒪_{K̂^{ur}}` の `(π^N)` に写る(降下の逆向き、易しい方)。 -/
theorem baseIntHom_mem_span_pow (K : PAdicLocalField p) (π : 𝒪[K.carrier]) (N : ℕ)
    {c : 𝒪[K.carrier]} (h : c ∈ Ideal.span {π ^ N}) :
    baseIntHom K c ∈ Ideal.span {uniformizerCompletionInt K π ^ N} := by
  rw [Ideal.mem_span_singleton'] at h ⊢
  obtain ⟨d, hd⟩ := h
  refine ⟨baseIntHom K d, ?_⟩
  have hBπ : baseIntHom K π = uniformizerCompletionInt K π := rfl
  rw [← hd, map_mul, map_pow, hBπ]

/-! ## 5. 主定理(環版) -/

set_option maxHeartbeats 1000000 in
/-- **★★★★★★★★★★★★★★★★★★★★★★(Λ6b)固定元は `𝒪_K` から来る**——
`σb = b` なら `b = a` なる `a ∈ 𝒪_K` がある。

`§3` で `b - a_N ∈ (π^N)` なる `a_N ∈ 𝒪_K` を各 `N` で取り、`§4` の降下で
`a_N` が `𝒪_K` の中で Cauchy であることを見る。`𝒪_K` は `𝔪` 進完備
(`isAdicComplete_valuationSubring` に `𝒪_K` の DVR 性を与えるだけ)なので
極限 `a` が取れ、`b - a ∈ (π^N)`(∀`N`)を `𝒪_{K̂^{ur}}` の `IsHausdorff` で
`0` にする。

★`set_option maxHeartbeats 1000000` は**配管の値段**である(部分型の上の
`map_sub` / `Submodule` の unification が既定を超える)。数学ではない。

退化の自己検査:`hσ` を落とすと**偽**(`σ = 1` は `hσ` を満たさず、
固定環は `𝒪_{K̂^{ur}}` 全体になる)。 -/
theorem exists_baseIntHom_eq_of_unramGalCompletionInt_eq (K : PAdicLocalField p)
    {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    {b : ↥(unramifiedCompletionInt K)} (hb : unramGalCompletionInt K σ b = b) :
    ∃ a : 𝒪[K.carrier], baseIntHom K a = b := by
  classical
  haveI hDVR : IsDiscreteValuationRing
      ((Valued.v : Valuation K.carrier NNReal).valuationSubring) :=
    isDiscreteValuationRing_carrierIntegers K
  haveI hACK : IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier] :=
    isAdicComplete_valuationSubring
  haveI hAC : IsAdicComplete (IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K))
      ↥(unramifiedCompletionInt K) := isAdicComplete_unramifiedCompletionInt K
  obtain ⟨π, hπ0, hπmax⟩ := exists_uniformizer_carrierIntegers K
  have hpowK : ∀ N : ℕ, IsLocalRing.maximalIdeal 𝒪[K.carrier] ^ N = Ideal.span {π ^ N} := by
    intro N; rw [hπmax, Ideal.span_singleton_pow]
  have hpowC := maximalIdeal_pow_unramifiedCompletionInt K hπ0 hπmax
  have hsub : ∀ {M N : ℕ}, M ≤ N →
      Ideal.span {uniformizerCompletionInt K π ^ N}
        ≤ Ideal.span {uniformizerCompletionInt K π ^ M} :=
    fun h => Ideal.span_singleton_le_span_singleton.mpr (pow_dvd_pow _ h)
  have key : ∀ N : ℕ, ∃ a : 𝒪[K.carrier],
      b - baseIntHom K a ∈ Ideal.span {uniformizerCompletionInt K π ^ N} := by
    intro N
    obtain ⟨a, b', -, hab⟩ :=
      exists_sub_baseIntHom_eq_uniformizer_pow_mul K hσ hπ0 hπmax N b hb
    exact ⟨a, by rw [hab, Ideal.mem_span_singleton']; exact ⟨b', mul_comm _ _⟩⟩
  choose A hA using key
  have hcauchy : ∀ {m n : ℕ}, m ≤ n →
      A m ≡ A n [SMOD (IsLocalRing.maximalIdeal 𝒪[K.carrier]) ^ m
        • (⊤ : Submodule 𝒪[K.carrier] 𝒪[K.carrier])] := by
    intro m n hmn
    rw [sModEq_pow_iff, hpowK]
    refine mem_span_pow_of_baseIntHom_mem K hπ0 m ?_
    have hrw : baseIntHom K (A m - A n)
        = (b - baseIntHom K (A n)) - (b - baseIntHom K (A m)) := by
      rw [map_sub]; ring
    rw [hrw]
    exact Submodule.sub_mem _ (hsub hmn (hA n)) (hA m)
  obtain ⟨a, ha⟩ := IsPrecomplete.prec inferInstance hcauchy
  refine ⟨a, ?_⟩
  have hzero : b - baseIntHom K a = 0 := by
    refine IsHausdorff.haus (I := IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K))
      inferInstance _ (fun N => ?_)
    rw [sModEq_pow_iff, sub_zero, hpowC]
    have h1 : A N - a ∈ Ideal.span {π ^ N} := by
      have h := (sModEq_pow_iff _ N _ _).mp (ha N)
      rwa [hpowK] at h
    have h3 : b - baseIntHom K a = (b - baseIntHom K (A N)) + baseIntHom K (A N - a) := by
      rw [map_sub]; ring
    rw [h3]
    exact Submodule.add_mem _ (hA N) (baseIntHom_mem_span_pow K π N h1)
  exact (sub_eq_zero.mp hzero).symm

/-- **★★★★★★★★★★★★★★★★★★(Λ6b)固定環の記述(同値の形)**——
`σb = b ⟺ b ∈ 𝒪_K の像`。 -/
theorem unramGalCompletionInt_eq_self_iff (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    (b : ↥(unramifiedCompletionInt K)) :
    unramGalCompletionInt K σ b = b ↔ ∃ a : 𝒪[K.carrier], baseIntHom K a = b := by
  refine ⟨fun hb => exists_baseIntHom_eq_of_unramGalCompletionInt_eq K hσ hb, ?_⟩
  rintro ⟨a, rfl⟩
  exact unramGalCompletionInt_baseIntHom K σ a

/-- 固定環の記述(集合の等式の形)。 -/
theorem setOf_unramGalCompletionInt_eq_self (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) :
    {b : ↥(unramifiedCompletionInt K) | unramGalCompletionInt K σ b = b}
      = ((baseIntHom K).range : Set ↥(unramifiedCompletionInt K)) := by
  ext b
  simpa using unramGalCompletionInt_eq_self_iff K hσ b

/-- 固定環の記述(`Subring` の等式の形)。`RingHom.eqLocus σ id` が固定環である。 -/
theorem eqLocus_unramGalCompletionInt (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) :
    RingHom.eqLocus (unramGalCompletionInt K σ : ↥(unramifiedCompletionInt K) →+*
        ↥(unramifiedCompletionInt K)) (RingHom.id _)
      = (baseIntHom K).range := by
  ext b
  simpa using unramGalCompletionInt_eq_self_iff K hσ b

/-- `𝒪_K` は固定環と**同型**である(`baseIntHom` が単射だから)。 -/
noncomputable def baseIntRangeEquiv (K : PAdicLocalField p) :
    𝒪[K.carrier] ≃+* (baseIntHom K).range :=
  RingEquiv.ofBijective (baseIntHom K).rangeRestrict
    ⟨fun _ _ h => baseIntHom_injective K (congrArg Subtype.val h),
      (baseIntHom K).rangeRestrict_surjective⟩

/-! ## 6. 単元版 -/

/-- **`𝒪_K^× → (𝒪_{K̂^{ur}})^×`**。 -/
noncomputable def baseUnitsHom (K : PAdicLocalField p) :
    (𝒪[K.carrier])ˣ →* (↥(unramifiedCompletionInt K))ˣ :=
  Units.map (baseIntHom K).toMonoidHom

@[simp] theorem baseUnitsHom_val (K : PAdicLocalField p) (v : (𝒪[K.carrier])ˣ) :
    ((baseUnitsHom K v : (↥(unramifiedCompletionInt K))ˣ) : ↥(unramifiedCompletionInt K))
      = baseIntHom K (v : 𝒪[K.carrier]) := rfl

/-- **★★★★★★★★★★★★★★(Λ6b)単元版の固定群**——
`σξ = ξ ⟺ ξ ∈ 𝒪_K^× の像`。

環版に落としてから、`baseIntHom` が**局所準同型**であること
(`isLocalHom_baseIntHom`)で「`a` は `𝒪_K` の単元」を回収する。

退化の自己検査:左辺を「固定される `0` でない元」に緩めると**偽**
(`π ≠ 0` は固定されるが `𝒪_K^×` の像には入らない)。 -/
theorem unramGalCompletionUnits_eq_self_iff (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    (ξ : (↥(unramifiedCompletionInt K))ˣ) :
    unramGalCompletionUnits K σ ξ = ξ ↔ ∃ v : (𝒪[K.carrier])ˣ, baseUnitsHom K v = ξ := by
  constructor
  · intro hξ
    have hval : unramGalCompletionInt K σ (ξ : ↥(unramifiedCompletionInt K))
        = (ξ : ↥(unramifiedCompletionInt K)) := by
      rw [← unramGalCompletionUnits_val K σ ξ, hξ]
    obtain ⟨a, ha⟩ := exists_baseIntHom_eq_of_unramGalCompletionInt_eq K hσ hval
    have hua : IsUnit a := by
      refine (isLocalHom_baseIntHom K).1 a ?_
      rw [ha]
      exact ξ.isUnit
    refine ⟨hua.unit, Units.ext ?_⟩
    rw [baseUnitsHom_val, IsUnit.unit_spec, ha]
  · rintro ⟨v, rfl⟩
    refine Units.ext ?_
    rw [unramGalCompletionUnits_val, baseUnitsHom_val, unramGalCompletionInt_baseIntHom]

/-- 単元版の固定群(集合の等式の形)。 -/
theorem setOf_unramGalCompletionUnits_eq_self (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) :
    {ξ : (↥(unramifiedCompletionInt K))ˣ | unramGalCompletionUnits K σ ξ = ξ}
      = ((baseUnitsHom K).range : Set (↥(unramifiedCompletionInt K))ˣ) := by
  ext ξ
  simpa using unramGalCompletionUnits_eq_self_iff K hσ ξ

/-! ## 7. 算術 Frobenius 版と主定理 -/

theorem arithFrobenius_eq_self_iff (K : PAdicLocalField p)
    (b : ↥(unramifiedCompletionInt K)) :
    unramGalCompletionInt K (arithFrobenius K) b = b
      ↔ ∃ a : 𝒪[K.carrier], baseIntHom K a = b :=
  unramGalCompletionInt_eq_self_iff K (arithFrobenius_residue K) b

theorem arithFrobenius_units_eq_self_iff (K : PAdicLocalField p)
    (ξ : (↥(unramifiedCompletionInt K))ˣ) :
    unramGalCompletionUnits K (arithFrobenius K) ξ = ξ
      ↔ ∃ v : (𝒪[K.carrier])ˣ, baseUnitsHom K v = ξ :=
  unramGalCompletionUnits_eq_self_iff K (arithFrobenius_residue K) ξ

/-- **★★★★★★★★★★★★★★★★★★★★★★★★(Λ6 Step 2 の入口)Dwork の補題(核)**——
`Gal(K^{ur}/K)` には、

* `𝓀_{K^{ur}}` 上で `z ↦ z^q` として作用し(= 算術 Frobenius)、
* **その 1 つの `σ` について**、`𝒪_{K̂^{ur}}` の中の `σ` の固定環が
  ちょうど `𝒪_K` の像である

元 `σ` がある。

★**量化の順が主張の中身**である。`σ` は `∃` の内側に閉じ込めてあり、
結論に自由なパラメータは出ない。

**証明**:`arithFrobenius K`(`ArithFrobeniusTopGen`)に
`setOf_unramGalCompletionInt_eq_self` を当てる。

退化の自己検査。

* **`σ` を任意の `unramGal K` に替えると偽**——`σ = 1` なら固定環は全体。
* `σ` を「位相的生成元」に弱めると**この証明は回らない**。
  ★主張自体が偽になるかは**未確認**。 -/
theorem exists_arithFrobenius_dworkFixedRing (K : PAdicLocalField p) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      {b : ↥(unramifiedCompletionInt K) | unramGalCompletionInt K σ b = b}
        = ((baseIntHom K).range : Set ↥(unramifiedCompletionInt K)) :=
  ⟨arithFrobenius K, arithFrobenius_residue K,
    setOf_unramGalCompletionInt_eq_self K (arithFrobenius_residue K)⟩

/-- **★★★★★★★★★★★★★★★★★★★★Dwork の補題(核・単元版)**——
`(𝒪_{K̂^{ur}})^×` の中の `σ` の固定群はちょうど `𝒪_K^×` の像である。 -/
theorem exists_arithFrobenius_dworkFixedUnits (K : PAdicLocalField p) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      {ξ : (↥(unramifiedCompletionInt K))ˣ | unramGalCompletionUnits K σ ξ = ξ}
        = ((baseUnitsHom K).range : Set (↥(unramifiedCompletionInt K))ˣ) :=
  ⟨arithFrobenius K, arithFrobenius_residue K,
    setOf_unramGalCompletionUnits_eq_self K (arithFrobenius_residue K)⟩

/-- **★★★★★★★★★★★★★★★★★★★★★★★★★★(Λ6)LEMMA 3.11 の全体**——
**同じ 1 つの `σ`** について、

1. `σ` は剰余体上 `z ↦ z^q`(算術 Frobenius)、
2. `σ` は `Gal(K^{ur}/K) ≅ Ẑ` の位相的生成元、
3. `b ↦ σb - b` は `𝒪_{K̂^{ur}}` 上**全射**(加法版、`DworkAdditive`)、
4. `ξ ↦ σξ/ξ` は `(𝒪_{K̂^{ur}})^×` 上**全射**(乗法版、`DworkMultiplicative`)、
5. 環の**核**:`{b | σb = b} = 𝒪_K の像`、
6. 単元の**核**:`{ξ | σξ = ξ} = 𝒪_K^× の像`

が同時に成り立つ。★6 つを別々の `∃` に分けたら無意味になる
(`σ` ごとに選び直せてしまう)。

★第 1059(加法版・乗法版)の逸脱記録「核は形式化していない」は、
本定理の 5・6 で**閉じた**。 -/
theorem exists_arithFrobenius_dworkKernelAndSurjective (K : PAdicLocalField p) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      (∀ N : ℕ, N ≠ 0 → σ ∈ unramLevelGeneratorSet K N) ∧
      (∀ c : ↥(unramifiedCompletionInt K), ∃ b : ↥(unramifiedCompletionInt K),
        unramGalCompletionInt K σ b - b = c) ∧
      (∀ u : (↥(unramifiedCompletionInt K))ˣ, ∃ ξ : (↥(unramifiedCompletionInt K))ˣ,
        unramGalCompletionUnits K σ ξ = ξ * u) ∧
      {b : ↥(unramifiedCompletionInt K) | unramGalCompletionInt K σ b = b}
        = ((baseIntHom K).range : Set ↥(unramifiedCompletionInt K)) ∧
      {ξ : (↥(unramifiedCompletionInt K))ˣ | unramGalCompletionUnits K σ ξ = ξ}
        = ((baseUnitsHom K).range : Set (↥(unramifiedCompletionInt K))ˣ) :=
  ⟨arithFrobenius K, arithFrobenius_residue K, arithFrobenius_mem_unramLevelGeneratorSet K,
    fun c => exists_arithFrobenius_sub_self_eq K c,
    fun u => exists_unramGalCompletionUnits_eq_mul K (arithFrobenius_residue K) u,
    setOf_unramGalCompletionInt_eq_self K (arithFrobenius_residue K),
    setOf_unramGalCompletionUnits_eq_self K (arithFrobenius_residue K)⟩

end ABC3.Found.PGC
