import ABC3.Found.PGC.DworkAdditive

/-!
# Dwork の補題(乗法版)—— `σξ/ξ = u` は `(𝒪_{K̂^{ur}})ˣ` で常に解ける(`sorry` 無し)

経路 Λ の節点 **Λ6 の 2 本目**。加法版(`DworkAdditive.lean`)の相方であり、
`Art_π` の `π` 非依存性が実際に使うのはこちらである:

```
∃ σ ∈ Gal(K^ur/K), (σ は剰余体で z ↦ z^q) ∧ ∀ u ∈ (𝒪_{K̂^ur})ˣ, ∃ ξ ∈ (𝒪_{K̂^ur})ˣ, σξ = ξ·u
```

★**量化の順が主張の中身**である。`∃ σ, ∀ u, ∃ ξ` であって `∀ u, ∃ σ, ∃ ξ` ではない
(後者は `u` ごとに Frobenius を選び直せるので空虚になる)。`σ` は `∃` の内側に
閉じ込めてあり、結論に自由なパラメータは出ない。

## ★原典がここを 1 行で畳んでいる

原典は Milne, *Class Field Theory* の LEMMA 3.11:

> The homomorphisms
> b ↦ b^σ − b : B → B;
> b ↦ b^σ/b : B^× → B^×,
> are surjective with kernels A and A^× respectively.

その証明は加法版(`B → B`)だけを書き、乗法版を

> The proof for B^× is similar.

の 1 行で畳む(`ResearchPaper/0_Source/Milne - Class Field Theory.txt` の 2501 行。
pdftotext は上付きを落とすので、その行は「The proof for B is similar.」と出る)。
畳まれているのは**乗法版 Dwork 丸ごと**である。測定器(`tools/hedge-index.mjs`)の
語彙が `similarly` しか見ておらず `is similar` を落としていたため、
2026-09-06 まで節点として数えられていなかった。

## 到達点

| | 主張 |
|---|---|
| 1 | `isUnit_of_sub_one_mem_maximalIdeal`:主単数は単元(局所環の一般論) |
| 2 | `exists_dworkMulStep_zero`:第 0 段(剰余体で `t̄^{q-1} = ū` を解く) |
| 3 | `exists_dworkMulStep_succ`:第 `n+1` 段(**加法版の可解性をそのまま使う**) |
| 4 | `exists_dworkMulStep`:2 つを 1 本の(仮定つき結論の)形にまとめる |
| 5 | `dworkPartialProd_spec`:部分積の不変式 `σ(ξ_N)·u_N = ξ_N·u` |
| 6 | `dworkMulResidual_sub_one_mem`:残差が `1 + 𝔪^n` を降りる |
| 7 | `exists_units_unramGalCompletionInt_eq_mul`:**`σξ = ξ·u` の可解性**(`σ` は仮引数) |
| 8 | **`exists_arithFrobenius_dworkMultiplicative`**:主定理(`σ` を `∃` の内側に) |
| 9 | `exists_arithFrobenius_isCoherent_dworkMultiplicative`:位相的生成元・加法版と**同じ 1 つの `σ`** |

## ★見立てとの差(2026-09-06 の実測)

段取り係の見立ては 3 つあった。実測は以下のとおり。

* 「**加法版の 1 段補題がそのまま効く**」——★**当たっている。しかも見立てより強い**。
  段取りは「`σe/e ≡ 1 + π^n(σa − a) mod 𝔪^{n+1}` を出して 1 段補題(`e^q − e = c̄`)を
  当てる」と読んでいたが、実際には**剰余体に降りる必要がなかった**。
  加法版の**主定理**(`exists_unramGalCompletionInt_sub_self_eq`、`σb − b = a` を
  誤差なしで解く)をそのまま呼ぶと、`t := 1 + π^{n+1} b` に対して

      σt = t + π^{n+1} a,   t·u = σt + π^{2(n+1)} a b

  という**恒等式**が出る(`ring` 1 行)。だから 1 段で `𝔪^{n+1} → 𝔪^{2(n+1)}` まで進む
  ——収束は 1 次ではなく **2 次**である。`n ≥ 1` では `2(n+1) ≥ n+2` なので
  必要な分より余る。★段取りが想定した「剰余体で `σa − a = −ā` を解く」段は**要らなかった**。
* 「**単数群の filtration `1 + 𝔪^n` が要る**」(加法版 docstring の見立て)——★**当たっている**。
  イデアルの減少列ではなく「`u_n − 1 ∈ (π^n)`」という主単数の条件で回っている。
* 「段 6(`ξ` が単元)が要る」——★**要った**。`σξ = ξ·u` だけからは `ξ = 0` を排除できない
  (`0` は解になってしまう)。`ξ ≡ ξ_1 mod 𝔪` と `ξ_1 = t_0` が単元であることから出す。
  ただし**安い**:4 行。

★**思ったより高かった所**:`𝒪_{K̂^{ur}}`(部分型の部分型)の上での `map_add` /
`map_mul` の unification が重く、`exists_dworkMulStep_succ` は既定の 200000 heartbeats を
超える。`set_option maxHeartbeats 1000000` を付けている(数学ではなく配管の値段)。

★**思ったより安かった所**:単元群と環の往復は、段取りが警告したほどには絡まなかった。
「先に `𝒪` の中で全部やって最後に単元化する」ではなく、**`Oˣ` の中で漸化式を回し、
不等式(イデアル所属)だけ `𝒪` に落とす**という分け方にしたら、`Units.ext` は
1 回しか要らなかった。逆元は `Oˣ` の中にしか現れない。

## 筋(6 段)

1. **第 0 段**。`u` は単元なので `ū ≠ 0`。`𝓀_{K̂^{ur}}` は代数閉なので
   `t̄^{q−1} = ū` が解け(`exists_pow_eq_completion`)、`t̄ ≠ 0` だから `t` は単元。
   `residue (σt) = t̄^q` より `t·u − σt ∈ 𝔪`。
2. **第 `n+1` 段**。`u − 1 = π^{n+1} a` に対し加法版で `σb − b = a` を解き、
   `t := 1 + π^{n+1} b` と置くと `t·u − σt = π^{2(n+1)} a b ∈ (π^{n+2})`。
   `t` は主単数なので単元。
3. 1 と 2 を `∃ t, (仮定 → 結論)` の**一様な**形にまとめる(`exists_dworkMulStep`)。
   ★こうしないと `choose` が仮定を引数に取る関数を作ってしまい、漸化式が組めない。
4. `choose T hT` で `T : ℕ → Oˣ → Oˣ` を取り、残差 `u_0 := u`,
   `u_{n+1} := (σt_n)⁻¹(t_n u_n)` と部分積 `ξ_N := ∏_{n<N} t_n` を定義する。
   不変式は `σ(ξ_N)·u_N = ξ_N·u`(帰納法 1 本、`linear_combination` で閉じる)。
5. `u_N − 1 ∈ (π^N)` と `t_N − 1 ∈ (π^N)` から `ξ_M − ξ_N ∈ (π^N)`。
   `IsPrecomplete.prec` で極限 `ξ` を取る。
6. `σξ − ξu = σ(ξ − ξ_N) − (ξ − ξ_N)u + (σξ_N − ξ_N u)` の 3 項がすべて `(π^N)` に
   入るので `IsHausdorff.haus` で `0`。`ξ ≡ ξ_1 mod 𝔪` から `ξ` は単元。

★ノルム・距離・収束は**一度も現れない**。`𝔪` 進完備性だけで閉じる(加法版と同じ)。

## 材料(すべて本プロジェクトの在庫)

| 要るもの | 出どころ |
|---|---|
| 算術 Frobenius | `ArithFrobeniusTopGen.arithFrobenius` / `arithFrobenius_residue` |
| `σb − b = a` の可解性 | `DworkAdditive.exists_unramGalCompletionInt_sub_self_eq` |
| `𝔪` 進完備性 / `𝔪^N = (π^N)` | `DworkAdditive.isAdicComplete_unramifiedCompletionInt` / `maximalIdeal_pow_unramifiedCompletionInt` |
| `σπ = π` / `(π^N)` の不変性 | `DworkAdditive.unramGalCompletionInt_uniformizerCompletionInt` / `unramGalCompletionInt_mem_span_pow` |
| 剰余体での `σ` の作用 | `DworkAdditive.residue_unramGalCompletionInt` |
| `z^n = a` の可解性(代数閉) | `UnramifiedResidueField.exists_pow_eq_completion` |
| `SModEq` の往復 | `DworkAdditive.sModEq_pow_iff` |

★**`frobeniusCompletionInt` は使っていない**。あれは位相的生成元一般から作られており、
剰余体上の作用は `z ↦ z^{q^k}` でよい。Dwork は**算術 Frobenius**でなければ回らない
(`UnramifiedResidueField.lean` の訂正欄)。

## ★設計上の注意(守ったこと)

* **既存の `Found/PGC/*.lean` を書き換えていない**(`DworkAdditive.lean` を import するだけ)。
* **`Skeleton` / `Interface` を触っていない**。`inertia` を経由せず、
  `Interface` の `SubgroupCorrespondence` / `ResidueCardinality` は主張にも証明にも
  現れない(Corollary 1.3 との循環を避けるため)。
* **`Abelianization` を使っていない**。
* **結論に自由なパラメータを出していない**——主定理の型は `K` にしか依存せず、
  `σ` も一意化元 `π` も `∃` / 証明の内側にある。
* **`sorry` 本体の `def` を作っていない**(D21)。

## 退化の自己検査

* **`σ` を任意の `unramGal K` に替えると偽**。`σ = 1` なら `σξ/ξ = 1` なので
  `u ≠ 1` は解けない。★だから `σ` は `∃` の内側に無ければならない。
* **`u` を `𝒪_{K̂^{ur}}` の一般の元(非単元を含む)に広げると偽**。`u = 0` は
  `σξ = ξ·0 = 0` を要求し、`ξ = 0` は単元でない(そして `ξ` が単元でない解を許すと
  下流の `Art_π` の議論が壊れる)。
* **剰余体が代数閉であることを落とすと第 0 段の `q−1` 乗根が取れない**。
  ★**主張自体が偽になるかは未確認 ——「偽」とは書かない**。
* **`σ` を「位相的生成元」に弱めるとこの証明は回らない**。第 0 段の剰余方程式が
  `t̄^{q^k−1} = ū`(`k ∈ Ẑ`)になり多項式方程式でなくなる。
  ★主張自体が偽になるかは**未確認**。
* `K^ur` を `K.closure` に替えると `𝔪` 進完備性が**偽**(`𝒪_{ℂ_p}` は DVR でない)。

## 逸脱(記録)

* **原典 LEMMA 3.11 の「kernel は `A^×`」の部分は形式化していない**。本ファイルが
  出しているのは**全射性だけ**である(`σξ/ξ = u` の可解性)。核の記述
  (`σξ = ξ ⟺ ξ ∈ 𝒪_K^×`)は下流(`Art_π` の `π` 非依存性)が使う形ではないので
  節点を立てていない。必要になったら立てること。
* **原典の証明路(逆極限と蛇の補題)を取っていない**。原典は `R/n^n` の完全列を
  `n` について帰納法で作り、逆極限(A.7 / A.8)へ渡す。本ファイルは加法版と同じく
  `IsAdicComplete` を位相の側から取り、逐次近似をイデアルの言葉だけで回した。
  ★Λ5 の逸脱(ノルム完備化と `𝔪` 進完備化の一致は未証明)を**必要としない**形である。
* 収束は原典の 1 次(各段で `n → n+1`)ではなく **2 次**(`n → 2n`)になっている。
  結論は同じなので下流への影響は無い。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found.GaloisRep
open scoped NNReal Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 0. 単元群への `σ` の作用と、主単数が単元であること -/

/-- `σ` の `𝒪_{K̂^{ur}}` の**単元群**への作用。逆元はここにしか現れない
(環 `𝒪_{K̂^{ur}}` の側では逆元を書かずに済ませる)。 -/
noncomputable def unramGalCompletionUnits (K : PAdicLocalField p) (σ : unramGal K) :
    (↥(unramifiedCompletionInt K))ˣ ≃* (↥(unramifiedCompletionInt K))ˣ :=
  Units.mapEquiv (unramGalCompletionInt K σ).toMulEquiv

@[simp] theorem unramGalCompletionUnits_val (K : PAdicLocalField p) (σ : unramGal K)
    (g : (↥(unramifiedCompletionInt K))ˣ) :
    ((unramGalCompletionUnits K σ g : (↥(unramifiedCompletionInt K))ˣ)
        : ↥(unramifiedCompletionInt K))
      = unramGalCompletionInt K σ (g : ↥(unramifiedCompletionInt K)) := rfl

theorem isUnit_unramGalCompletionInt_val (K : PAdicLocalField p) (σ : unramGal K)
    (g : (↥(unramifiedCompletionInt K))ˣ) :
    IsUnit (unramGalCompletionInt K σ (g : ↥(unramifiedCompletionInt K))) :=
  ⟨unramGalCompletionUnits K σ g, rfl⟩

/-- **主単数は単元**——`x - 1 ∈ 𝔪` なら `x` は単元(局所環の一般論)。

★乗法版が加法版と違うのはここから先で、収束を測るのがイデアルの減少列 `(π^n)` ではなく
**主単数の filtration `1 + 𝔪^n`** になる。その出入りに毎回要る。 -/
theorem isUnit_of_sub_one_mem_maximalIdeal {R : Type*} [CommRing R] [IsLocalRing R] {x : R}
    (h : x - 1 ∈ IsLocalRing.maximalIdeal R) : IsUnit x := by
  by_contra hx
  have hxm : x ∈ IsLocalRing.maximalIdeal R := (IsLocalRing.mem_maximalIdeal x).mpr hx
  have hone : (1 : R) ∈ IsLocalRing.maximalIdeal R := by
    have h2 := Submodule.sub_mem _ hxm h
    rwa [sub_sub_cancel] at h2
  exact (IsLocalRing.maximalIdeal.isMaximal R).ne_top ((Ideal.eq_top_iff_one _).mpr hone)

/-! ## 1. 第 0 段 —— 剰余体で `t̄^{q-1} = ū` を解く

加法版の第 1 段は Artin–Schreier(`ē^q - ē = c̄`)だったが、乗法版は
`𝓀_{K̂^{ur}}` が代数閉であることから `q-1` 乗根を取る。 -/

/-- **★★★★★★★★★★(Λ6)乗法版の第 0 段**——単元 `u` に対し、単元 `t` で
`t·u - σt ∈ 𝔪` なるものが取れる。

`ū ≠ 0` なので `t̄^{q-1} = ū` が解け(`𝓀_{K̂^{ur}}` は代数閉)、`t̄ ≠ 0` だから
持ち上げた `t` は単元。`residue (σt) = t̄^q` なので
`residue (t·u - σt) = t̄·t̄^{q-1} - t̄^q = 0`。

退化の自己検査。

* `u` が単元でないと `ū = 0` になり、`t̄` も `0` になるので `t` が単元にならない。
* `hσ`(`σ` が剰余体で `q` 乗)を落とすと**偽**——`σ = 1` なら `t·u - t` は
  一般に `𝔪` に入らない。 -/
theorem exists_dworkMulStep_zero (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    (w : (↥(unramifiedCompletionInt K))ˣ) :
    ∃ t : (↥(unramifiedCompletionInt K))ˣ,
      (t : ↥(unramifiedCompletionInt K)) * (w : ↥(unramifiedCompletionInt K))
          - unramGalCompletionInt K σ (t : ↥(unramifiedCompletionInt K))
        ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K) := by
  have hq1 : 1 < Nat.card 𝓀[K.carrier] := one_lt_residueCard K
  have hwne : IsLocalRing.residue _ (w : ↥(unramifiedCompletionInt K)) ≠ 0 := by
    rw [Ne, IsLocalRing.residue_eq_zero_iff, IsLocalRing.mem_maximalIdeal, mem_nonunits_iff]
    exact not_not.mpr w.isUnit
  obtain ⟨s, hs⟩ := exists_pow_eq_completion K
    (IsLocalRing.residue _ (w : ↥(unramifiedCompletionInt K)))
    (n := Nat.card 𝓀[K.carrier] - 1) (by omega)
  obtain ⟨t0, rfl⟩ := IsLocalRing.residue_surjective s
  have ht0ne : IsLocalRing.residue _ t0 ≠ 0 := by
    intro h0
    rw [h0, zero_pow (by omega)] at hs
    exact hwne hs.symm
  have ht0 : IsUnit t0 := by
    by_contra hcon
    exact ht0ne ((IsLocalRing.residue_eq_zero_iff t0).mpr
      ((IsLocalRing.mem_maximalIdeal t0).mpr hcon))
  refine ⟨ht0.unit, ?_⟩
  rw [IsUnit.unit_spec, ← IsLocalRing.residue_eq_zero_iff, map_sub, map_mul,
    residue_unramGalCompletionInt K hσ t0, ← hs, ← pow_succ']
  rw [show Nat.card 𝓀[K.carrier] - 1 + 1 = Nat.card 𝓀[K.carrier] by omega, sub_self]

/-! ## 2. 第 `n+1` 段 —— 加法版の可解性がそのまま効く

★段取りの見立て(「剰余体に降りて `σa - a = -ā` を解く」)より**強いことが言える**。
加法版の主定理は `σb - b = a` を**誤差なしで**解くので、`t := 1 + π^{n+1} b` に対し

    σt = t + π^{n+1} a,    t·u = σt + π^{2(n+1)} a b

が**恒等式**として出る(`ring` 1 行)。1 段で `𝔪^{n+1} → 𝔪^{2(n+1)}`、すなわち
収束は 2 次である。 -/

set_option maxHeartbeats 1000000 in
/-- **★★★★★★★★★★★★(Λ6)乗法版の第 `n+1` 段**——`u - 1 ∈ (π^{n+1})` なる単元 `u` に
対し、単元 `t` で `t - 1 ∈ (π^{n+1})` かつ `t·u - σt ∈ (π^{n+2})` なるものが取れる。

`u - 1 = a·π^{n+1}` と書き、**加法版**で `σb - b = a` を解いて `t := 1 + π^{n+1} b`。
`σπ = π` から `σt = 1 + π^{n+1}(b + a)` で、差は `π^{2(n+1)} a b`。

退化の自己検査。

* `hw`(`u ≡ 1 mod π^{n+1}`)を落とすと `t` が主単数にならず、`t·u - σt` の評価が崩れる。
* `hπmax` を落とすと `π` が自由なパラメータになり、`t` の単元性(`π` が非単元)が言えない。
* `σπ = π` を落とすと `σt = 1 + (σπ)^{n+1}(b+a)` になり恒等式が壊れる。 -/
theorem exists_dworkMulStep_succ (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    {π : 𝒪[K.carrier]} (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (n : ℕ) (w : (↥(unramifiedCompletionInt K))ˣ)
    (hw : (w : ↥(unramifiedCompletionInt K)) - 1
      ∈ Ideal.span {uniformizerCompletionInt K π ^ (n + 1)}) :
    ∃ t : (↥(unramifiedCompletionInt K))ˣ,
      ((t : ↥(unramifiedCompletionInt K)) - 1
        ∈ Ideal.span {uniformizerCompletionInt K π ^ (n + 1)}) ∧
      ((t : ↥(unramifiedCompletionInt K)) * (w : ↥(unramifiedCompletionInt K))
          - unramGalCompletionInt K σ (t : ↥(unramifiedCompletionInt K))
        ∈ Ideal.span {uniformizerCompletionInt K π ^ (n + 2)}) := by
  obtain ⟨a, ha⟩ := Ideal.mem_span_singleton'.mp hw
  obtain ⟨b, hb⟩ := exists_unramGalCompletionInt_sub_self_eq K hσ a
  have hσb : unramGalCompletionInt K σ b = b + a := by rw [← hb]; ring
  have hσt0 : unramGalCompletionInt K σ (1 + uniformizerCompletionInt K π ^ (n + 1) * b)
      = 1 + uniformizerCompletionInt K π ^ (n + 1) * (b + a) := by
    rw [map_add, map_one, map_mul, map_pow,
      unramGalCompletionInt_uniformizerCompletionInt K σ π, hσb]
  have ht0mem : (1 + uniformizerCompletionInt K π ^ (n + 1) * b) - 1
      ∈ Ideal.span {uniformizerCompletionInt K π ^ (n + 1)} := by
    rw [Ideal.mem_span_singleton']
    exact ⟨b, by ring⟩
  have ht0unit : IsUnit (1 + uniformizerCompletionInt K π ^ (n + 1) * b) := by
    refine isUnit_of_sub_one_mem_maximalIdeal ?_
    rw [maximalIdeal_unramifiedCompletionInt_eq_span K hπ0 hπmax, Ideal.mem_span_singleton']
    exact ⟨uniformizerCompletionInt K π ^ n * b, by ring⟩
  have hwv : (w : ↥(unramifiedCompletionInt K))
      = 1 + a * uniformizerCompletionInt K π ^ (n + 1) := by linear_combination -ha
  refine ⟨ht0unit.unit, by rw [IsUnit.unit_spec]; exact ht0mem, ?_⟩
  rw [IsUnit.unit_spec, hσt0, hwv, Ideal.mem_span_singleton']
  exact ⟨uniformizerCompletionInt K π ^ n * a * b, by ring⟩

/-! ## 3. 2 つの段を 1 本にまとめる

★`choose` で漸化式の材料を取り出すには、存在文の**仮定が `t` に依存しない**形、
すなわち `∃ t, (仮定 → 結論)` でなければならない。仮定を `∀` の外に出すと
`choose` が「証明を引数に取る関数」を作ってしまい、`ℕ → Oˣ → Oˣ` にならない。 -/

set_option maxHeartbeats 1000000 in
/-- **★★★★★★★★★★★★★★(Λ6)乗法版の 1 段(一様形)**——`n = 0` は剰余体の
`q-1` 乗根(§1)、`n ≥ 1` は加法版(§2)。

`n = 0` では仮定 `u - 1 ∈ (π^0) = (1)` が自明に成り立ち、結論の第 1 項
`t - 1 ∈ (π^0)` も自明。だから 2 つの段が同じ型に収まる。 -/
theorem exists_dworkMulStep (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    {π : 𝒪[K.carrier]} (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (n : ℕ) (w : (↥(unramifiedCompletionInt K))ˣ) :
    ∃ t : (↥(unramifiedCompletionInt K))ˣ,
      (w : ↥(unramifiedCompletionInt K)) - 1
          ∈ Ideal.span {uniformizerCompletionInt K π ^ n} →
        ((t : ↥(unramifiedCompletionInt K)) - 1
            ∈ Ideal.span {uniformizerCompletionInt K π ^ n} ∧
          (t : ↥(unramifiedCompletionInt K)) * (w : ↥(unramifiedCompletionInt K))
              - unramGalCompletionInt K σ (t : ↥(unramifiedCompletionInt K))
            ∈ Ideal.span {uniformizerCompletionInt K π ^ (n + 1)}) := by
  classical
  cases n with
  | zero =>
    obtain ⟨t, ht⟩ := exists_dworkMulStep_zero K hσ w
    refine ⟨t, fun _ => ⟨?_, ?_⟩⟩
    · rw [pow_zero, Ideal.span_singleton_one]; exact Submodule.mem_top
    · rw [zero_add, pow_one, ← maximalIdeal_unramifiedCompletionInt_eq_span K hπ0 hπmax]
      exact ht
  | succ m =>
    by_cases hw : (w : ↥(unramifiedCompletionInt K)) - 1
        ∈ Ideal.span {uniformizerCompletionInt K π ^ (m + 1)}
    · obtain ⟨t, h1, h2⟩ := exists_dworkMulStep_succ K hσ hπ0 hπmax m w hw
      exact ⟨t, fun _ => ⟨h1, h2⟩⟩
    · exact ⟨1, fun h => absurd h hw⟩

/-! ## 4. 逐次近似 —— 残差と部分積

★加法版の `dworkPartialSum`(部分**和**)に対応するのが `dworkPartialProd`(部分**積**)。
残差 `u_n` は加法版の `C^[n] c` に当たるが、乗法版では逆元が要るので `Oˣ` の中で回す。 -/

/-- 逐次近似の残差 `u_n`:`u_0 = u`、`u_{n+1} = (σt_n)⁻¹·(t_n·u_n)`。 -/
noncomputable def dworkMulResidual (K : PAdicLocalField p) (σ : unramGal K)
    (T : ℕ → (↥(unramifiedCompletionInt K))ˣ → (↥(unramifiedCompletionInt K))ˣ)
    (u : (↥(unramifiedCompletionInt K))ˣ) : ℕ → (↥(unramifiedCompletionInt K))ˣ
  | 0 => u
  | n + 1 =>
      (unramGalCompletionUnits K σ (T n (dworkMulResidual K σ T u n)))⁻¹
        * (T n (dworkMulResidual K σ T u n) * dworkMulResidual K σ T u n)

@[simp] theorem dworkMulResidual_zero (K : PAdicLocalField p) (σ : unramGal K)
    (T : ℕ → (↥(unramifiedCompletionInt K))ˣ → (↥(unramifiedCompletionInt K))ˣ)
    (u : (↥(unramifiedCompletionInt K))ˣ) : dworkMulResidual K σ T u 0 = u := rfl

@[simp] theorem dworkMulResidual_succ (K : PAdicLocalField p) (σ : unramGal K)
    (T : ℕ → (↥(unramifiedCompletionInt K))ˣ → (↥(unramifiedCompletionInt K))ˣ)
    (u : (↥(unramifiedCompletionInt K))ˣ) (n : ℕ) :
    dworkMulResidual K σ T u (n + 1)
      = (unramGalCompletionUnits K σ (T n (dworkMulResidual K σ T u n)))⁻¹
          * (T n (dworkMulResidual K σ T u n) * dworkMulResidual K σ T u n) := rfl

/-- 逐次近似の**部分積** `ξ_N = ∏_{n<N} t_n`。★加法版の部分和に対応する所。 -/
noncomputable def dworkPartialProd (K : PAdicLocalField p) (σ : unramGal K)
    (T : ℕ → (↥(unramifiedCompletionInt K))ˣ → (↥(unramifiedCompletionInt K))ˣ)
    (u : (↥(unramifiedCompletionInt K))ˣ) : ℕ → (↥(unramifiedCompletionInt K))ˣ
  | 0 => 1
  | N + 1 => dworkPartialProd K σ T u N * T N (dworkMulResidual K σ T u N)

@[simp] theorem dworkPartialProd_zero (K : PAdicLocalField p) (σ : unramGal K)
    (T : ℕ → (↥(unramifiedCompletionInt K))ˣ → (↥(unramifiedCompletionInt K))ˣ)
    (u : (↥(unramifiedCompletionInt K))ˣ) : dworkPartialProd K σ T u 0 = 1 := rfl

@[simp] theorem dworkPartialProd_succ (K : PAdicLocalField p) (σ : unramGal K)
    (T : ℕ → (↥(unramifiedCompletionInt K))ˣ → (↥(unramifiedCompletionInt K))ˣ)
    (u : (↥(unramifiedCompletionInt K))ˣ) (N : ℕ) :
    dworkPartialProd K σ T u (N + 1)
      = dworkPartialProd K σ T u N * T N (dworkMulResidual K σ T u N) := rfl

/-- 残差の定義そのもの:`σ(t_n)·u_{n+1} = t_n·u_n`(単元群の中)。 -/
theorem unramGalCompletionUnits_mul_dworkMulResidual_succ (K : PAdicLocalField p)
    (σ : unramGal K)
    (T : ℕ → (↥(unramifiedCompletionInt K))ˣ → (↥(unramifiedCompletionInt K))ˣ)
    (u : (↥(unramifiedCompletionInt K))ˣ) (n : ℕ) :
    unramGalCompletionUnits K σ (T n (dworkMulResidual K σ T u n))
        * dworkMulResidual K σ T u (n + 1)
      = T n (dworkMulResidual K σ T u n) * dworkMulResidual K σ T u n := by
  rw [dworkMulResidual_succ, mul_inv_cancel_left]

/-- 同じ関係を環 `𝒪_{K̂^{ur}}` の言葉に落としたもの。★以降、逆元は出てこない。 -/
theorem unramGalCompletionInt_mul_dworkMulResidual_succ (K : PAdicLocalField p)
    (σ : unramGal K)
    (T : ℕ → (↥(unramifiedCompletionInt K))ˣ → (↥(unramifiedCompletionInt K))ˣ)
    (u : (↥(unramifiedCompletionInt K))ˣ) (n : ℕ) :
    unramGalCompletionInt K σ ((T n (dworkMulResidual K σ T u n) : ↥(unramifiedCompletionInt K)))
        * (dworkMulResidual K σ T u (n + 1) : ↥(unramifiedCompletionInt K))
      = (T n (dworkMulResidual K σ T u n) : ↥(unramifiedCompletionInt K))
          * (dworkMulResidual K σ T u n : ↥(unramifiedCompletionInt K)) := by
  calc unramGalCompletionInt K σ
        ((T n (dworkMulResidual K σ T u n) : ↥(unramifiedCompletionInt K)))
          * (dworkMulResidual K σ T u (n + 1) : ↥(unramifiedCompletionInt K))
      = ((unramGalCompletionUnits K σ (T n (dworkMulResidual K σ T u n))
            * dworkMulResidual K σ T u (n + 1) : (↥(unramifiedCompletionInt K))ˣ)
          : ↥(unramifiedCompletionInt K)) := by
        rw [Units.val_mul, unramGalCompletionUnits_val]
    _ = ((T n (dworkMulResidual K σ T u n) * dworkMulResidual K σ T u n
          : (↥(unramifiedCompletionInt K))ˣ) : ↥(unramifiedCompletionInt K)) := by
        rw [unramGalCompletionUnits_mul_dworkMulResidual_succ]
    _ = _ := Units.val_mul _ _

set_option maxHeartbeats 1000000 in
/-- **★★★★★★★★★★★★(Λ6)逐次近似の不変式(乗法版)**:`σ(ξ_N)·u_N = ξ_N·u`。

`N` の帰納法 1 本。`σ` を積に分配するのは `map_mul` と `Units.val_mul` だけで、
あとは `linear_combination σ(ξ_N)·hstep + t_N·ih` で閉じる。

★加法版の `dworkPartialSum_spec` に対応。あちらが `σπ = π` を要としたのに対し、
こちらの要は**残差の定義**(`σ(t_N)·u_{N+1} = t_N·u_N`)である。 -/
theorem dworkPartialProd_spec (K : PAdicLocalField p) (σ : unramGal K)
    (T : ℕ → (↥(unramifiedCompletionInt K))ˣ → (↥(unramifiedCompletionInt K))ˣ)
    (u : (↥(unramifiedCompletionInt K))ˣ) (N : ℕ) :
    unramGalCompletionInt K σ (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K))
        * (dworkMulResidual K σ T u N : ↥(unramifiedCompletionInt K))
      = (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K))
          * (u : ↥(unramifiedCompletionInt K)) := by
  induction N with
  | zero => rw [dworkPartialProd_zero, dworkMulResidual_zero, Units.val_one, map_one]
  | succ N ih =>
    have hstep := unramGalCompletionInt_mul_dworkMulResidual_succ K σ T u N
    rw [dworkPartialProd_succ, Units.val_mul, map_mul]
    linear_combination
      unramGalCompletionInt K σ (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K))
          * hstep
        + (T N (dworkMulResidual K σ T u N) : ↥(unramifiedCompletionInt K)) * ih

set_option maxHeartbeats 1000000 in
/-- **★★★★★★★★★★(Λ6)残差は主単数の filtration を降りる**:`u_n - 1 ∈ (π^n)`。

★これが乗法版の収束の管理である(加法版は `(π^N)` というイデアルの減少列だった)。
`σ(t_n)` が単元であることを使って `σ(t_n)·(u_{n+1} - 1) = t_n·u_n - σ(t_n)` の
所属から `u_{n+1} - 1` の所属へ移る(`Ideal.unit_mul_mem_iff_mem`)。 -/
theorem dworkMulResidual_sub_one_mem (K : PAdicLocalField p) (σ : unramGal K)
    {π : 𝒪[K.carrier]}
    {T : ℕ → (↥(unramifiedCompletionInt K))ˣ → (↥(unramifiedCompletionInt K))ˣ}
    (hT : ∀ (n : ℕ) (w : (↥(unramifiedCompletionInt K))ˣ),
      (w : ↥(unramifiedCompletionInt K)) - 1
          ∈ Ideal.span {uniformizerCompletionInt K π ^ n} →
        ((T n w : ↥(unramifiedCompletionInt K)) - 1
            ∈ Ideal.span {uniformizerCompletionInt K π ^ n} ∧
          (T n w : ↥(unramifiedCompletionInt K)) * (w : ↥(unramifiedCompletionInt K))
              - unramGalCompletionInt K σ (T n w : ↥(unramifiedCompletionInt K))
            ∈ Ideal.span {uniformizerCompletionInt K π ^ (n + 1)}))
    (u : (↥(unramifiedCompletionInt K))ˣ) (n : ℕ) :
    (dworkMulResidual K σ T u n : ↥(unramifiedCompletionInt K)) - 1
      ∈ Ideal.span {uniformizerCompletionInt K π ^ n} := by
  induction n with
  | zero =>
    rw [pow_zero, Ideal.span_singleton_one]
    exact Submodule.mem_top
  | succ n ih =>
    obtain ⟨-, h2⟩ := hT n (dworkMulResidual K σ T u n) ih
    refine (Ideal.unit_mul_mem_iff_mem _
      (isUnit_unramGalCompletionInt_val K σ (T n (dworkMulResidual K σ T u n)))).mp ?_
    have heq : unramGalCompletionInt K σ
          (T n (dworkMulResidual K σ T u n) : ↥(unramifiedCompletionInt K))
          * ((dworkMulResidual K σ T u (n + 1) : ↥(unramifiedCompletionInt K)) - 1)
        = (T n (dworkMulResidual K σ T u n) : ↥(unramifiedCompletionInt K))
            * (dworkMulResidual K σ T u n : ↥(unramifiedCompletionInt K))
          - unramGalCompletionInt K σ
              (T n (dworkMulResidual K σ T u n) : ↥(unramifiedCompletionInt K)) := by
      rw [mul_sub, mul_one, unramGalCompletionInt_mul_dworkMulResidual_succ]
    rw [heq]
    exact h2

set_option maxHeartbeats 1000000 in
/-- **部分積は `(π^N)` を法として安定**:`ξ_{N+k} - ξ_N ∈ (π^N)`。

★`IsPrecomplete.prec` が要求する Cauchy 条件はこれ(と差の符号の入れ替え)である。
加法版の `dworkPartialSum_sub_mem` に対応するが、差の分解が
`ξ_{N+k+1} - ξ_N = ξ_{N+k}(t_{N+k} - 1) + (ξ_{N+k} - ξ_N)` になる点が違う。
★`N = 0` の段(`t_0 - 1` は小さくない)を落とさないため、条件は
`t_n - 1 ∈ (π^n)`(`n = 0` では自明)という形でなければならない。 -/
theorem dworkPartialProd_sub_mem (K : PAdicLocalField p) (σ : unramGal K)
    {π : 𝒪[K.carrier]}
    {T : ℕ → (↥(unramifiedCompletionInt K))ˣ → (↥(unramifiedCompletionInt K))ˣ}
    (hT : ∀ (n : ℕ) (w : (↥(unramifiedCompletionInt K))ˣ),
      (w : ↥(unramifiedCompletionInt K)) - 1
          ∈ Ideal.span {uniformizerCompletionInt K π ^ n} →
        ((T n w : ↥(unramifiedCompletionInt K)) - 1
            ∈ Ideal.span {uniformizerCompletionInt K π ^ n} ∧
          (T n w : ↥(unramifiedCompletionInt K)) * (w : ↥(unramifiedCompletionInt K))
              - unramGalCompletionInt K σ (T n w : ↥(unramifiedCompletionInt K))
            ∈ Ideal.span {uniformizerCompletionInt K π ^ (n + 1)}))
    (u : (↥(unramifiedCompletionInt K))ˣ) (N k : ℕ) :
    (dworkPartialProd K σ T u (N + k) : ↥(unramifiedCompletionInt K))
        - (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K))
      ∈ Ideal.span {uniformizerCompletionInt K π ^ N} := by
  induction k with
  | zero => simp
  | succ k ih =>
    obtain ⟨c, hc⟩ := Ideal.mem_span_singleton'.mp
      (hT (N + k) (dworkMulResidual K σ T u (N + k))
        (dworkMulResidual_sub_one_mem K σ hT u (N + k))).1
    have hfirst : (dworkPartialProd K σ T u (N + k) : ↥(unramifiedCompletionInt K))
        * ((T (N + k) (dworkMulResidual K σ T u (N + k)) : ↥(unramifiedCompletionInt K)) - 1)
        ∈ Ideal.span {uniformizerCompletionInt K π ^ N} := by
      rw [Ideal.mem_span_singleton']
      exact ⟨(dworkPartialProd K σ T u (N + k) : ↥(unramifiedCompletionInt K)) * c
        * uniformizerCompletionInt K π ^ k, by rw [← hc, pow_add]; ring⟩
    have heq : (dworkPartialProd K σ T u (N + (k + 1)) : ↥(unramifiedCompletionInt K))
          - (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K))
        = (dworkPartialProd K σ T u (N + k) : ↥(unramifiedCompletionInt K))
            * ((T (N + k) (dworkMulResidual K σ T u (N + k))
                : ↥(unramifiedCompletionInt K)) - 1)
          + ((dworkPartialProd K σ T u (N + k) : ↥(unramifiedCompletionInt K))
              - (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K))) := by
      rw [← Nat.add_assoc, dworkPartialProd_succ, Units.val_mul]; ring
    rw [heq]
    exact Submodule.add_mem _ hfirst ih

/-! ## 5. 可解性 -/

set_option maxHeartbeats 1000000 in
/-- **★★★★★★★★★★★★★★★★★★(Λ6)`σξ = ξ·u` は `(𝒪_{K̂^{ur}})ˣ` で解ける**
(`σ` を仮引数に取った形)。

`§3` を `choose` して `T` を取り、残差 `u_N` と部分積 `ξ_N`(`§4`)を作る。
`dworkPartialProd_sub_mem` から `IsPrecomplete.prec` で極限 `ξ` を取り、

```
σξ - ξu = σ(ξ - ξ_N) - (ξ - ξ_N)u + (σξ_N - ξ_N u)
```

の 3 項がすべて `(π^N)` に入ることを見て `IsHausdorff.haus` で `0` にする。
第 3 項は不変式 `σ(ξ_N)·u_N = ξ_N·u` から `σ(ξ_N)·(1 - u_N)` に等しく、
`u_N - 1 ∈ (π^N)`(`§4`)で小さい。

★**`ξ` が単元であること**は自動ではない(`ξ = 0` は `σξ = ξu` を満たす)。
`ξ ≡ ξ_1 mod 𝔪` と `ξ_1 = t_0` が単元であることから出す。

★ノルム・距離・収束は使わない。`𝔪` 進完備性だけで閉じる。

退化の自己検査:`hσ` を落とすと**偽**(`σ = 1` は `hσ` を満たさず、
`σξ = ξ` なので `u ≠ 1` は解けない)。 -/
theorem exists_units_unramGalCompletionInt_eq_mul (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    (u : (↥(unramifiedCompletionInt K))ˣ) :
    ∃ ξ : (↥(unramifiedCompletionInt K))ˣ,
      unramGalCompletionInt K σ (ξ : ↥(unramifiedCompletionInt K))
        = (ξ : ↥(unramifiedCompletionInt K)) * (u : ↥(unramifiedCompletionInt K)) := by
  classical
  haveI hAC : IsAdicComplete (IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K))
      ↥(unramifiedCompletionInt K) := isAdicComplete_unramifiedCompletionInt K
  obtain ⟨π, hπ0, hπmax⟩ := exists_uniformizer_carrierIntegers K
  choose T hT using exists_dworkMulStep K hσ hπ0 hπmax
  have hpow := maximalIdeal_pow_unramifiedCompletionInt K hπ0 hπmax
  have hcauchy : ∀ {m n : ℕ}, m ≤ n →
      (dworkPartialProd K σ T u m : ↥(unramifiedCompletionInt K))
        ≡ (dworkPartialProd K σ T u n : ↥(unramifiedCompletionInt K))
        [SMOD (IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K)) ^ m
          • (⊤ : Submodule ↥(unramifiedCompletionInt K) ↥(unramifiedCompletionInt K))] := by
    intro m n hmn
    rw [sModEq_pow_iff, hpow]
    have h := dworkPartialProd_sub_mem K σ hT u m (n - m)
    rw [Nat.add_sub_cancel' hmn] at h
    simpa using Submodule.neg_mem _ h
  obtain ⟨x, hx⟩ := IsPrecomplete.prec inferInstance hcauchy
  have hxunit : IsUnit x := by
    have h1 : x - (dworkPartialProd K σ T u 1 : ↥(unramifiedCompletionInt K))
        ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K) := by
      have h := (sModEq_pow_iff _ 1 _ _).mp (hx 1)
      rw [pow_one] at h
      simpa using Submodule.neg_mem _ h
    by_contra hcon
    have hξ1 : (dworkPartialProd K σ T u 1 : ↥(unramifiedCompletionInt K))
        ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K) := by
      have h2 := Submodule.sub_mem _ ((IsLocalRing.mem_maximalIdeal x).mpr hcon) h1
      rwa [sub_sub_cancel] at h2
    exact (mem_nonunits_iff.mp ((IsLocalRing.mem_maximalIdeal _).mp hξ1))
      (dworkPartialProd K σ T u 1).isUnit
  refine ⟨hxunit.unit, ?_⟩
  rw [IsUnit.unit_spec]
  refine sub_eq_zero.mp ?_
  refine IsHausdorff.haus (I := IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K))
    inferInstance _ (fun N => ?_)
  rw [sModEq_pow_iff, sub_zero, hpow]
  have h1 : x - (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K))
      ∈ Ideal.span {uniformizerCompletionInt K π ^ N} := by
    have h := (sModEq_pow_iff _ N _ _).mp (hx N)
    rw [hpow] at h
    simpa using Submodule.neg_mem _ h
  have h3 := unramGalCompletionInt_mem_span_pow K σ π N h1
  have h4 : (x - (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K)))
      * (u : ↥(unramifiedCompletionInt K))
      ∈ Ideal.span {uniformizerCompletionInt K π ^ N} := Ideal.mul_mem_right _ _ h1
  obtain ⟨d, hd⟩ := Ideal.mem_span_singleton'.mp
    (dworkMulResidual_sub_one_mem K σ hT u N)
  have h2 : unramGalCompletionInt K σ
        (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K))
      - (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K))
          * (u : ↥(unramifiedCompletionInt K))
      ∈ Ideal.span {uniformizerCompletionInt K π ^ N} := by
    rw [Ideal.mem_span_singleton']
    refine ⟨-(unramGalCompletionInt K σ
      (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K)) * d), ?_⟩
    have hspec := dworkPartialProd_spec K σ T u N
    linear_combination -hspec - unramGalCompletionInt K σ
      (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K)) * hd
  have hdecomp : unramGalCompletionInt K σ x - x * (u : ↥(unramifiedCompletionInt K))
      = unramGalCompletionInt K σ
            (x - (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K)))
          - (x - (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K)))
              * (u : ↥(unramifiedCompletionInt K))
        + (unramGalCompletionInt K σ
              (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K))
            - (dworkPartialProd K σ T u N : ↥(unramifiedCompletionInt K))
                * (u : ↥(unramifiedCompletionInt K))) := by
    rw [map_sub]; ring
  rw [hdecomp]
  exact Submodule.add_mem _ (Submodule.sub_mem _ h3 h4) h2

set_option maxHeartbeats 1000000 in
/-- 単元群の中で書いた形:`(Units.map σ) ξ = ξ · u`。 -/
theorem exists_unramGalCompletionUnits_eq_mul (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    (u : (↥(unramifiedCompletionInt K))ˣ) :
    ∃ ξ : (↥(unramifiedCompletionInt K))ˣ, unramGalCompletionUnits K σ ξ = ξ * u := by
  obtain ⟨ξ, hξ⟩ := exists_units_unramGalCompletionInt_eq_mul K hσ u
  refine ⟨ξ, Units.ext ?_⟩
  rw [unramGalCompletionUnits_val, Units.val_mul]
  exact hξ

set_option maxHeartbeats 1000000 in
/-- 環の言葉(`IsUnit`)で書いた形。下流が `Oˣ` を持ちたくないとき用。 -/
theorem exists_isUnit_unramGalCompletionInt_eq_mul (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    {u : ↥(unramifiedCompletionInt K)} (hu : IsUnit u) :
    ∃ ξ : ↥(unramifiedCompletionInt K), IsUnit ξ ∧
      unramGalCompletionInt K σ ξ = ξ * u := by
  obtain ⟨ξ, hξ⟩ := exists_units_unramGalCompletionInt_eq_mul K hσ hu.unit
  exact ⟨(ξ : ↥(unramifiedCompletionInt K)), ξ.isUnit, by rw [hξ, IsUnit.unit_spec]⟩

set_option maxHeartbeats 1000000 in
/-- **原典 LEMMA 3.11 の言い方**:`ξ ↦ ξ^σ/ξ` は `(𝒪_{K̂^{ur}})ˣ` 上**全射**。

★原典が主張している「kernel は `A^×`」の部分は形式化していない
(モジュール docstring の逸脱記録)。 -/
theorem surjective_unramGalCompletionUnits_div_self (K : PAdicLocalField p) :
    Function.Surjective (fun ξ : (↥(unramifiedCompletionInt K))ˣ =>
      unramGalCompletionUnits K (arithFrobenius K) ξ * ξ⁻¹) := by
  intro u
  obtain ⟨ξ, hξ⟩ := exists_unramGalCompletionUnits_eq_mul K (arithFrobenius_residue K) u
  refine ⟨ξ, ?_⟩
  show unramGalCompletionUnits K (arithFrobenius K) ξ * ξ⁻¹ = u
  rw [hξ, mul_comm (ξ : (↥(unramifiedCompletionInt K))ˣ) u, mul_assoc,
    mul_inv_cancel, mul_one]

/-! ## 6. 主定理 -/

set_option maxHeartbeats 1000000 in
/-- **★★★★★★★★★★★★★★★★★★★★★★★★(Λ6)Dwork の補題(乗法版)**——
`Gal(K^ur/K)` には、

* `𝓀_{K^{ur}}` 上で `z ↦ z^q` として作用し(= 算術 Frobenius)、
* **その 1 つの `σ` について**、任意の単元 `u ∈ (𝒪_{K̂^ur})ˣ` に対し
  `σξ = ξ·u` なる単元 `ξ ∈ (𝒪_{K̂^ur})ˣ` が存在する

元 `σ` がある。★原典(Milne, CFT の LEMMA 3.11)が
「The proof for B^× is similar.」の 1 行で畳んでいた所である。

★**量化の順が主張の中身**である。`∃ σ, ∀ u, ∃ ξ` を主張しており、
`∀ u, ∃ σ, ∃ ξ` ではない(後者は `u` ごとに Frobenius を選び直せるので空虚)。

**証明**:`arithFrobenius K`(`ArithFrobeniusTopGen`)に
`exists_units_unramGalCompletionInt_eq_mul` を当てる。

退化の自己検査。

* **`σ` を任意の `unramGal K` に替えると偽**——`σ = 1` なら `σξ/ξ = 1`。
* **`u` を非単元まで広げると偽**——`u = 0` は `ξ = 0` を強いる。
* 剰余体が代数閉であることを落とすと第 0 段の `q-1` 乗根が取れない
  (★主張自体が偽になるかは**未確認**)。 -/
theorem exists_arithFrobenius_dworkMultiplicative (K : PAdicLocalField p) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      ∀ u : (↥(unramifiedCompletionInt K))ˣ, ∃ ξ : (↥(unramifiedCompletionInt K))ˣ,
        unramGalCompletionInt K σ (ξ : ↥(unramifiedCompletionInt K))
          = (ξ : ↥(unramifiedCompletionInt K)) * (u : ↥(unramifiedCompletionInt K)) :=
  ⟨arithFrobenius K, arithFrobenius_residue K,
    fun u => exists_units_unramGalCompletionInt_eq_mul K (arithFrobenius_residue K) u⟩

set_option maxHeartbeats 1000000 in
/-- **★★★★★★★★★★★★★★★★★★★★★★同じ `σ` が位相的生成元でも加法版の解でもある形**。

Λ5(`Ẑ` との同定)と Λ6 の加法版・乗法版が**同じ 1 つの元**を指せるように、
4 つの性質を 1 つの `∃` に並べた形も出しておく。★4 つを別々の `∃` に分けたら無意味になる。

これが `Art_π` の `π` 非依存性(Λ6 本体)が消費する形である。 -/
theorem exists_arithFrobenius_isCoherent_dworkMultiplicative (K : PAdicLocalField p) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      (∀ N : ℕ, N ≠ 0 → σ ∈ unramLevelGeneratorSet K N) ∧
      (∀ c : ↥(unramifiedCompletionInt K), ∃ b : ↥(unramifiedCompletionInt K),
        unramGalCompletionInt K σ b - b = c) ∧
      ∀ u : (↥(unramifiedCompletionInt K))ˣ, ∃ ξ : (↥(unramifiedCompletionInt K))ˣ,
        unramGalCompletionInt K σ (ξ : ↥(unramifiedCompletionInt K))
          = (ξ : ↥(unramifiedCompletionInt K)) * (u : ↥(unramifiedCompletionInt K)) :=
  ⟨arithFrobenius K, arithFrobenius_residue K, arithFrobenius_mem_unramLevelGeneratorSet K,
    fun c => exists_arithFrobenius_sub_self_eq K c,
    fun u => exists_units_unramGalCompletionInt_eq_mul K (arithFrobenius_residue K) u⟩

end ABC3.Found.PGC
