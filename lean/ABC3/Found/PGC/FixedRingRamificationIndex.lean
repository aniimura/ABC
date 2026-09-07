import ABC3.Found.PGC.QuotientRamIndexAverage

/-!
# 固定環の分岐指数 `e(K′/K′′) = |H|`(★原典が名前を付けていない基盤ノード)

★★**これは原典が独立の主張として立てていない基盤ノードである。**
Yoshida は §6.2 で `K′/K` を "totally ramified" と置き、`e(K′/K′′) = |H|` を
その中に畳んでいる(名前も番号も付いていない)。したがって本ファイルには
`.src`(`ABC3.Meta.Source`)を**置かない** —— 存在しない `sectionId` を書くことになるからである。

**なぜ立てたか**: Y9(`ABC3/Found/PGC/QuotientRamIndexAverage.lean`、Yoshida 2008 Lemma 6.8)が
唯一これだけを仮定に置いて着地したからである:

```
(hram : ∀ c : C, addVal B (algebraMap C B c) = (Nat.card H : ℕ∞) * addVal C c)
```

★本ファイルはこれを**証明する**(`addVal_map_eq_card_mul`)。
さらに Y9 の Lemma 6.8 から `hram` を外した形
`card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing` を最後に置いた
(★Y9 のファイルは書き換えていない。差し替えるかどうかは上流の判断に委ねる)。

## 設定

`B` = `O_{K′}`(DVR、`G` が作用)、`C` = `O_{K′′}` = `B^H`(DVR)、`H ≤ G` 有限。
`algebraMap C B` は単射で、`H` はその像を各点固定する(`hinvC`)。

* `hfix` : `B^H ⊆ algebraMap C B '' C` —— 「`C` は固定環である」の実質。
* `hres` : `∀ b : B, ∃ c : C, b − algebraMap C B c ∈ m_B` —— ★**完全分岐 (`f = 1`)**。
  剰余体が伸びないこと。★これを落とすと主張は**偽**(不分岐なら `e = 1`、`f = |H|`)。
* `hadj` : `Algebra.adjoin C {π} = ⊤` —— `O_{K′} = O_{K′′}[π′]`(原典 Lemma 5.11 の帰結)。
* `[FaithfulSMul G B]` —— ★これを落とすと主張は**偽**(`H` が自明に作用すれば `C = B`、`e = 1`)。

## 証明の骨(道 A + Nakayama)

**§1(抽象核・付値だけ)** `addVal_map_eq_mul_addVal` : 離散付値の延長は素元での値だけで決まる。
`c = u ϖ^k` と分解して `addVal B` を当てるだけ。★分岐も Galois も群も出てこない。

**§3(道 A の前半)** `ϖ₀ := ∏_{τ∈H} τ•π`(ノルム)は `H` 不変だから `C` に入り、
`addVal B ϖ₀ = |H|`(付値は群作用で不変、Y1 の `addVal_smul` + Y9 の `addVal_prod`)。
§1 と合わせて `|H| = e · v_C(ϖ₀)`、すなわち **`e ∣ |H|` と `e ≤ |H|`**。

**§4(残る山 `|H| ≤ e`)** ★本体の段取りは「`v_C(ϖ₀) = 1`(`ϖ₀` が `C` の素元)を示す」だったが、
それは直接には出ない。代わりに次の 2 段で閉じた。

* §4a `exists_pow_eq_sum`(★**抽象核 —— 群作用が出てこない**):
  完全分岐 (`hres`) と有限性 (`hfg`) から `π` は `C` 上**次数 `e` の monic 関係**
  `π^e = ∑_{i<e} c_i π^i` を満たす。
  段取り: `M := span_C{π^i : i < e}` と置き、「`v(b) ≥ j` なら `M` の元を引いて `v ≥ e` にできる」を
  `j` の**有限降下帰納**で示す(各段で `hres` が 1 桁分の係数を `C` から取ってくる)。
  これで `⊤ = M ⊔ m_C • ⊤` となり、**中山の補題**
  (`Submodule.le_of_le_smul_of_le_jacobson_bot`)で `M = ⊤`。
  ★`m_C • ⊤ = π^e B` は `algebraMap C B ϖ` と `π^e` の同伴(付値が等しい)から出る。
* §4b `card_le_of_pow_eq_sum`(具体層): 係数 `c_i` は `C` の元だから `H` で不変。
  よって**相異なる**共役 `τ•π`(`τ ∈ H`)がすべて同じ monic 関係の根になる。
  次数 `e` の多項式は `|H|` 個の相異なる根を持てないので **`|H| ≤ e`**
  (Y9 の `prod_X_sub_C_dvd_of_forall_isRoot` + `Polynomial.natDegree_le_of_dvd`)。

**§6** `e ≤ |H|` と `|H| ≤ e` から `e = |H|`、§1 に戻して主結論。

★★**道 B(Artin の `[B:C] = |H|` と `e·f = n`)は使わなかった。**
`Ideal.ramificationIdx` と `addVal` を繋ぐ橋も、`e·f = n` も要らない。
中山の補題だけで足りる(実測: §4a は 1 往復で通った)。

## ★★逸脱の記録

1. ★**`hres`(完全分岐)を「剰余体が伸びない」の形で仮定に置いた。**
   原典は §6.2 で `K′/K` が totally ramified だと述べており、そこから `K′/K′′` も
   completely ramified であることは(局所体の理論では)標準だが、
   ★**この木にも mathlib にも「`K′/K` 完全分岐 ⟹ `K′/K′′` 完全分岐」は無い**。
   本ファイルは `K′/K′′` の側で仮定した。★弱めても強めてもいない読み替えである。
2. ★**`hadj : Algebra.adjoin C {π} = ⊤`(`O_{K′} = O_{K′′}[π′]`)を仮定に置いた。**
   原典 Lemma 5.11 は `O_{K′} = O[π′]`(`O = O_K` 上)を与える。
   `O ⊆ O_{K′′}` だから前者は後者から従う —— その導出も
   `adjoin_eq_top_of_image` として本ファイルに入れてある(仮定は `O ⊆ O_{K′′}` の
   `hAC : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a` だけ)。
3. `C` を `B` の部分環として構成せず `[Algebra C B]` + `hinj` + `hfix` で受けた(Y9 と同じ)。
   ★`MulSemiringAction G C` は**要らなかった** —— 主結論が使うのは
   `hinvC : ∀ ρ ∈ H, ∀ x : C, ρ • algebraMap C B x = algebraMap C B x` の 1 本だけである
   (Y9 の `hcomp` + `hHtriv` から 1 行で出る。§7 で実際にそうしている)。
4. `H` の正規性も `G/H` の商群も使っていない(Y9 と同じ理由)。

## 退化の自己検査

* ★**`e = 0` を許すと §1 が偽になる**(`c = 0` のとき左辺 `⊤`、右辺 `0 * ⊤ = 0`)。
  `addVal_map_eq_mul_addVal` は `he : addVal B (algebraMap C B ϖ) ≠ 0` を**明示的に要求する**。
  ★これは `algebraMap` の単射性からは出ない(反例: `C = ℤ_p ↪ B = ℚ_p[[t]]` は単射だが
  `p` は `B` の単元)。主結論では §2 `isUnit_of_isUnit_map`(単元は固定環に降りる)で塞いだ。
* ★**`H` が無限だと `Nat.card H = 0` で右辺が潰れ、空虚に真**になる。`[Fintype ↥H]` で塞がっている
  (Y9 の `pos_natCard_subgroup` と同じ穴)。
* ★**完全分岐 `hres` を落とすと偽**(不分岐拡大 `W(F_q) ⊆ W(F_{q^n})` では `e = 1 ≠ n = |H|`)。
* ★**`[FaithfulSMul G B]` を落とすと偽**(`H` が `B` に自明に作用すれば `C = B`、`e = 1`)。
  §5 `smul_conj_injective` で「共役 `τ•π` が相異なる」に変換して使う。
* ★`algebraMap C B` の単射性を落とすと `hinj` を使う 2 箇所(§2 と `e ≠ ⊤`)が崩れる。

★`open Polynomial` しているので `Polynomial.C` と型変数 `C`(= `O_{K′′}`)が衝突する。
多項式の定数は**すべて `Polynomial.C` と修飾して**書いてある。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing
open Polynomial

/-! ## §1 抽象核 -/

/-- ★★★**抽象核** —— 離散付値の延長は素元での値だけで決まる。

`C ⊆ B` を離散付値環の(単射とは限らない)代数とし、`ϖ` を `C` の素元、
`e := addVal B (algebraMap C B ϖ)` と置くと、**すべての** `c : C` について

`addVal B (algebraMap C B c) = e * addVal C c`。

★**分岐・Galois・群作用の語彙は 1 つも出てこない。**
`c = u ϖ^k`(`u` は単元)と分解し、`algebraMap` が単元を単元に写すことを使うだけ。

★退化検査: `he : e ≠ 0` は**落とせない**。`c = 0` のとき左辺は `⊤`、右辺は `0 * ⊤ = 0` で偽になる。
★`algebraMap` の単射性からは `e ≠ 0` は出ない(`e ≠ ⊤` しか出ない)。 -/
theorem addVal_map_eq_mul_addVal {C B : Type*} [CommRing C] [IsDomain C]
    [IsDiscreteValuationRing C] [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    [Algebra C B] {ϖ : C} (hϖ : Irreducible ϖ) (he : addVal B (algebraMap C B ϖ) ≠ 0) (c : C) :
    addVal B (algebraMap C B c) = addVal B (algebraMap C B ϖ) * addVal C c := by
  rcases eq_or_ne c 0 with rfl | hc
  · simp only [map_zero, addVal_zero, ENat.mul_top he]
  obtain ⟨k, u, rfl⟩ := eq_unit_mul_pow_irreducible hc hϖ
  rw [addVal_def' u hϖ k, map_mul, map_pow, addVal_mul, addVal_pow,
    addVal_eq_zero_iff.mpr (u.isUnit.map (algebraMap C B)), zero_add, nsmul_eq_mul, mul_comm]

/-! ## §2 単元は降りる -/

/-- **単元は固定環に降りる** —— `algebraMap C B c` が `B` の単元なら `c` は `C` の単元。

`algebraMap C B c` は `H` 不変(`hinvC`)なので、その逆元も `H` 不変であり、
`hfix`(`B^H ⊆ C` の像)によって `C` から来る。あとは `hinj` で `c * c' = 1`。

★これが `e ≠ 0`(§1 の仮定)の供給源である。 -/
theorem isUnit_of_isUnit_map {C B : Type*} [CommRing C] [CommRing B] [Algebra C B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G}
    (hinj : Function.Injective (algebraMap C B))
    (hinvC : ∀ ρ ∈ H, ∀ c : C, ρ • algebraMap C B c = algebraMap C B c)
    (hfix : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    {c : C} (h : IsUnit (algebraMap C B c)) : IsUnit c := by
  obtain ⟨v, hv⟩ := h
  set y : B := ((v⁻¹ : Bˣ) : B) with hy
  have hcy : algebraMap C B c * y = 1 := by rw [← hv, hy, v.mul_inv]
  have hinvy : ∀ ρ ∈ H, ρ • y = y := by
    intro ρ hρ
    have h1 : algebraMap C B c * (ρ • y) = 1 := by
      rw [← hinvC ρ hρ c, ← smul_mul', hcy, smul_one]
    calc ρ • y = (ρ • y) * (algebraMap C B c * y) := by rw [hcy, mul_one]
      _ = (algebraMap C B c * (ρ • y)) * y := by ring
      _ = y := by rw [h1, one_mul]
  obtain ⟨c', hc'⟩ := hfix y hinvy
  refine IsUnit.of_mul_eq_one c' (hinj ?_)
  rw [map_mul, hc', hcy, map_one]

/-! ## §3 共役の積 -/

/-- **ノルム元 `ϖ₀ := ∏_{τ∈H} τ•π` は `H` 不変** —— `ρ ∈ H` は `H` を平行移動で置換する。
★Y9 の `map_conjProd_self`(多項式版)の元版。 -/
theorem smul_prod_conj {B : Type*} [CommRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    {H : Subgroup G} [Fintype H] (π : B) {ρ : G} (hρ : ρ ∈ H) :
    ρ • (∏ τ : H, (τ : G) • π) = ∏ τ : H, (τ : G) • π := by
  have h0 : ρ • (∏ τ : H, (τ : G) • π) = ∏ τ : H, ρ • ((τ : G) • π) :=
    map_prod (MulSemiringAction.toRingHom G B ρ) _ _
  calc ρ • (∏ τ : H, (τ : G) • π)
      = ∏ τ : H, ρ • ((τ : G) • π) := h0
    _ = ∏ τ : H, (((Equiv.mulLeft (⟨ρ, hρ⟩ : H)) τ : H) : G) • π :=
        Finset.prod_congr rfl fun τ _ => by rw [← mul_smul]; rfl
    _ = ∏ τ : H, (τ : G) • π :=
        Equiv.prod_comp (Equiv.mulLeft (⟨ρ, hρ⟩ : H)) (fun τ : H => (τ : G) • π)

/-- **`v_{K′}(∏_{τ∈H} τ•π) = |H|`** —— 付値は群作用で不変(`addVal_smul`)で、
素元の付値は `1`。積の付値は和(Y9 の `addVal_prod`)。 -/
theorem addVal_prod_conj {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G} [Fintype H] {π : B}
    (hπ : Irreducible π) :
    addVal B (∏ τ : H, (τ : G) • π) = (Nat.card H : ℕ∞) := by
  rw [addVal_prod]
  simp [addVal_smul, addVal_uniformizer hπ, Nat.card_eq_fintype_card]


/-! ## §4a 抽象核 —— 完全分岐なら素元は次数 `e` の monic 関係を満たす -/

/-- ★★★**抽象核(段 4a)** —— 完全分岐なら素元 `π` は `C` 上**次数 `e` の monic 関係**を満たす。

`e := addVal B (algebraMap C B ϖ)`(`ϖ ∈ m_C`)、`hres` が「剰余体が伸びない」、
`hfg` が「`B` は `C`-加群として有限生成」であるとき

`π^e = ∑_{i<e} algebraMap C B (c i) * π^i`。

★**群作用も Galois も出てこない。** 段取りは
1. `M := span_C {π^i : i < e}`。
2. `key` : `v(b) ≥ j` なら `M` の元を引いて `v ≥ e` にできる(`j` の有限降下帰納)。
   各段で `hres` が「1 桁分の係数」を `C` から取ってくる。
3. `⊤ = M ⊔ m_C • ⊤`(`m_C • ⊤ = π^e B` は `algebraMap C B ϖ ~ᵤ π^e` から)。
4. **中山の補題** `Submodule.le_of_le_smul_of_le_jacobson_bot` で `M = ⊤`。
5. `π^e ∈ M` の展開が求める関係。

★`e = 0` でも主張は意味を持つ(そのとき仮定は矛盾していて `1 = 0` が出る)。
別に場合分けする必要は無い。 -/
theorem exists_pow_eq_sum {C B : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C]
    [CommRing B] [IsDomain B] [IsDiscreteValuationRing B] [Algebra C B]
    {π : B} (hπ : Irreducible π) {ϖ : C} (hϖ : ϖ ∈ maximalIdeal C) {e : ℕ}
    (he : addVal B (algebraMap C B ϖ) = (e : ℕ∞))
    (hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hfg : (⊤ : Submodule C B).FG) :
    ∃ c : Fin e → C, π ^ e = ∑ i : Fin e, algebraMap C B (c i) * π ^ (i : ℕ) := by
  classical
  set M : Submodule C B := Submodule.span C (Set.range fun i : Fin e => π ^ (i : ℕ)) with hM
  have hmax : maximalIdeal B = Ideal.span {π} := (irreducible_iff_uniformizer π).mp hπ
  have hmem : ∀ b : B, b ∈ maximalIdeal B ↔ π ∣ b := by
    intro b; rw [hmax, Ideal.mem_span_singleton]
  have key : ∀ d j : ℕ, e ≤ j + d → ∀ b : B, π ^ j ∣ b → ∃ m ∈ M, π ^ e ∣ (b - m) := by
    intro d
    induction d with
    | zero =>
        intro j hj b hb
        exact ⟨0, Submodule.zero_mem _, by
          rw [sub_zero]; exact dvd_trans (pow_dvd_pow π (by omega)) hb⟩
    | succ d ih =>
        intro j hj b hb
        by_cases hje : e ≤ j
        · exact ⟨0, Submodule.zero_mem _, by
            rw [sub_zero]; exact dvd_trans (pow_dvd_pow π hje) hb⟩
        replace hje : j < e := Nat.not_le.mp hje
        obtain ⟨b', rfl⟩ := hb
        obtain ⟨c, hc⟩ := hres b'
        obtain ⟨w, hw⟩ := (hmem _).1 hc
        have hnext : π ^ (j + 1) ∣ π ^ j * b' - algebraMap C B c * π ^ j := by
          refine ⟨w, ?_⟩
          have hfac : π ^ j * b' - algebraMap C B c * π ^ j
              = π ^ j * (b' - algebraMap C B c) := by ring
          rw [hfac, hw, pow_succ]
          ring
        obtain ⟨m, hmM, hm⟩ := ih (j + 1) (by omega) _ hnext
        refine ⟨algebraMap C B c * π ^ j + m, Submodule.add_mem _ ?_ hmM, ?_⟩
        · rw [hM, ← Algebra.smul_def]
          exact Submodule.smul_mem _ _ (Submodule.subset_span ⟨⟨j, hje⟩, rfl⟩)
        · have hrw : π ^ j * b' - (algebraMap C B c * π ^ j + m)
              = (π ^ j * b' - algebraMap C B c * π ^ j) - m := by ring
          rw [hrw]
          exact hm
  have hassoc : Associated (algebraMap C B ϖ) (π ^ e) := by
    rw [← addVal_eq_iff_associated, he, addVal_pow, addVal_uniformizer hπ]
    simp
  have hjac : maximalIdeal C ≤ Ideal.jacobson (⊥ : Ideal C) :=
    le_of_eq (IsLocalRing.jacobson_eq_maximalIdeal ⊥ bot_ne_top).symm
  have htop : (⊤ : Submodule C B) ≤ M ⊔ (maximalIdeal C) • (⊤ : Submodule C B) := by
    intro b _
    obtain ⟨m, hmM, hdvd⟩ := key e 0 (by omega) b (by simp)
    obtain ⟨z, hz⟩ := hassoc.dvd.trans hdvd
    have hbm : b - m ∈ (maximalIdeal C) • (⊤ : Submodule C B) := by
      rw [hz, ← Algebra.smul_def]
      exact Submodule.smul_mem_smul hϖ Submodule.mem_top
    have hsum := Submodule.add_mem_sup hmM hbm
    rwa [add_sub_cancel] at hsum
  have hMtop : (⊤ : Submodule C B) ≤ M :=
    Submodule.le_of_le_smul_of_le_jacobson_bot hfg hjac htop
  have hpow : π ^ e ∈ M := hMtop Submodule.mem_top
  rw [hM, Submodule.mem_span_range_iff_exists_fun] at hpow
  obtain ⟨c, hc⟩ := hpow
  exact ⟨c, by rw [← hc]; exact Finset.sum_congr rfl fun i _ => Algebra.smul_def _ _⟩

/-! ## §4b 具体層 —— 相異なる共役が同じ次数 `e` の monic 関係を満たすので `|H| ≤ e` -/

/-- **段 4b(具体層)** —— 係数が `C` から来る次数 `e` の monic 関係を `π` が満たし、
共役 `τ•π`(`τ ∈ H`)が相異なるなら `|H| ≤ e`。

係数は `H` で不変(`hinvC`)なので、`τ` を関係式に当てると `τ•π` も同じ関係を満たす。
`g := X^e − ∑ C(c i) X^i` は monic で次数 `e`、その根が `|H|` 個あるので
`∏_{τ∈H}(X − τ•π) ∣ g`(Y9 の `prod_X_sub_C_dvd_of_forall_isRoot`)、次数を較べて `|H| ≤ e`。 -/
theorem card_le_of_pow_eq_sum {C B : Type*} [CommRing C] [CommRing B] [IsDomain B] [Algebra C B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G} [Fintype H] {π : B} {e : ℕ}
    {c : Fin e → C}
    (hinvC : ∀ ρ ∈ H, ∀ x : C, ρ • algebraMap C B x = algebraMap C B x)
    (hrel : π ^ e = ∑ i : Fin e, algebraMap C B (c i) * π ^ (i : ℕ))
    (hdist : ∀ τ₁ τ₂ : H, (τ₁ : G) • π = (τ₂ : G) • π → τ₁ = τ₂) :
    Nat.card H ≤ e := by
  classical
  set g : B[X] := X ^ e - ∑ i : Fin e, Polynomial.C (algebraMap C B (c i)) * X ^ (i : ℕ) with hg
  have hlt : (∑ i : Fin e, Polynomial.C (algebraMap C B (c i)) * X ^ (i : ℕ)).degree
      < ((X : B[X]) ^ e).degree := by
    rw [Polynomial.degree_X_pow]
    refine lt_of_le_of_lt (Polynomial.degree_sum_le _ _) ?_
    refine (Finset.sup_lt_iff (by exact WithBot.bot_lt_coe e)).2 fun i _ => ?_
    exact lt_of_le_of_lt (Polynomial.degree_C_mul_X_pow_le _ _)
      (by exact_mod_cast Nat.cast_lt.2 i.isLt)
  have hmonic : g.Monic := (Polynomial.monic_X_pow e).sub_of_left hlt
  have hdeg : g.natDegree = e := by
    refine Polynomial.natDegree_eq_of_degree_eq_some ?_
    rw [hg, Polynomial.degree_sub_eq_left_of_degree_lt hlt, Polynomial.degree_X_pow]
  have hroot : ∀ τ : H, g.eval ((τ : G) • π) = 0 := by
    intro τ
    have hsm : (τ : G) • (π ^ e)
        = ∑ i : Fin e, algebraMap C B (c i) * ((τ : G) • π) ^ (i : ℕ) := by
      rw [hrel, Finset.smul_sum]
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [smul_mul', smul_pow', hinvC (τ : G) τ.2]
    rw [smul_pow'] at hsm
    rw [hg]
    simp only [Polynomial.eval_sub, Polynomial.eval_pow, Polynomial.eval_X,
      Polynomial.eval_finsetSum, Polynomial.eval_mul, Polynomial.eval_C]
    rw [hsm, sub_self]
  have hdvd : (∏ τ : H, (X - Polynomial.C ((τ : G) • π))) ∣ g :=
    prod_X_sub_C_dvd_of_forall_isRoot (fun i _ j _ h => hdist i j h) (fun τ _ => hroot τ)
  have hnd : (∏ τ : H, (X - Polynomial.C ((τ : G) • π))).natDegree = Nat.card H := by
    rw [Polynomial.natDegree_prod_of_monic _ _ fun τ _ => Polynomial.monic_X_sub_C _]
    simp [Nat.card_eq_fintype_card]
  calc Nat.card H = (∏ τ : H, (X - Polynomial.C ((τ : G) • π))).natDegree := hnd.symm
    _ ≤ g.natDegree := Polynomial.natDegree_le_of_dvd hdvd hmonic.ne_zero
    _ = e := hdeg

/-! ## §5 配管 —— 共役の相異性と有限生成性 -/

/-- `τ` が `C` の像と `π` を固定し `B = C[π]` なら、`τ` は `B` 全体を固定する。 -/
theorem smul_eq_self_of_adjoin_eq_top {C B : Type*} [CommRing C] [CommRing B] [Algebra C B]
    {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hadj : Algebra.adjoin C ({π} : Set B) = ⊤) {τ : G}
    (hτC : ∀ x : C, τ • algebraMap C B x = algebraMap C B x) (hτπ : τ • π = π) (b : B) :
    τ • b = b := by
  have hb : b ∈ Algebra.adjoin C ({π} : Set B) := by rw [hadj]; trivial
  induction hb using Algebra.adjoin_induction with
  | mem x hx => rw [Set.mem_singleton_iff.1 hx]; exact hτπ
  | algebraMap x => exact hτC x
  | add x y _ _ hx hy => rw [smul_add, hx, hy]
  | mul x y _ _ hx hy => rw [smul_mul', hx, hy]

/-- **共役の相異性** —— `τ₁•π = τ₂•π` なら `τ₁ = τ₂`。
★`[FaithfulSMul G B]` を使う唯一の場所。これを落とすと主結論は偽になる
(`H` が自明に作用すれば `C = B` で `e = 1`)。 -/
theorem smul_conj_injective {C B : Type*} [CommRing C] [CommRing B] [Algebra C B]
    {G : Type*} [Group G] [MulSemiringAction G B] [FaithfulSMul G B] {H : Subgroup G} {π : B}
    (hadj : Algebra.adjoin C ({π} : Set B) = ⊤)
    (hinvC : ∀ ρ ∈ H, ∀ x : C, ρ • algebraMap C B x = algebraMap C B x)
    (τ₁ τ₂ : H) (h : (τ₁ : G) • π = (τ₂ : G) • π) : τ₁ = τ₂ := by
  have h1 : ((τ₂ : G)⁻¹ * (τ₁ : G)) • π = π := by rw [mul_smul, h, inv_smul_smul]
  have hmem : ((τ₂ : G)⁻¹ * (τ₁ : G)) ∈ H := H.mul_mem (H.inv_mem τ₂.2) τ₁.2
  have h2 : ((τ₂ : G)⁻¹ * (τ₁ : G)) = 1 :=
    eq_of_smul_eq_smul (α := B) fun x => by
      rw [smul_eq_self_of_adjoin_eq_top hadj (hinvC _ hmem) h1 x, one_smul]
  exact Subtype.ext (inv_mul_eq_one.mp h2).symm

/-- `π` は `C` 上整 —— `∏_{τ∈H}(X − τ•π)` は monic でその係数は `H` 不変(Y9 の
`smul_coeff_conjProd`)だから `C` から来る(`hfix`)。`Polynomial.lifts_and_degree_eq_and_monic`
で `C[X]` に持ち上げる。★§6 で `hfg`(`B` が `C`-加群として有限生成)を作るために要る。 -/
theorem isIntegral_of_fixed {C B : Type*} [CommRing C] [CommRing B] [Nontrivial B] [Algebra C B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G} [Fintype H] {π : B}
    (hfix : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b) :
    IsIntegral C π := by
  classical
  have hmonicf : (∏ τ : H, (X - Polynomial.C ((τ : G) • π))).Monic :=
    Polynomial.monic_prod_of_monic _ _ fun τ _ => Polynomial.monic_X_sub_C _
  have hlifts : (∏ τ : H, (X - Polynomial.C ((τ : G) • π))) ∈
      Polynomial.lifts (algebraMap C B) := by
    rw [Polynomial.lifts_iff_coeff_lifts]
    exact fun n => hfix _ fun ρ hρ => smul_coeff_conjProd ρ hρ n
  obtain ⟨q, hq, -, hqm⟩ := Polynomial.lifts_and_degree_eq_and_monic hlifts hmonicf
  refine ⟨q, hqm, ?_⟩
  rw [Polynomial.eval₂_eq_eval_map, hq, eval_prod_X_sub_C]
  exact Finset.prod_eq_zero (Finset.mem_univ (1 : H)) (by simp)

/-! ## §6 主結論 -/

/-- ★★★★**主結論 —— `e(K′/K′′) = |H|` の掛け算形**。

`∀ c : C, addVal B (algebraMap C B c) = |H| * addVal C c`。

★**これが Y9 `card_mul_ramIndex_eq_sum_ramIndex` の仮定 `hram` そのものである。**
Y9 は「導いていない。仮定である」と書いた。ここで証明した。

仮定の役割:
* `hπ` : `π` は `B` の素元(共役の積がノルム元になる)。
* `hinj` : `algebraMap C B` 単射(`e ≠ ⊤` に要る)。
* `hinvC` : `H` は `C` の像を各点固定する(`C ⊆ B^H`)。
* `hfix` : `B^H ⊆ C` の像(`C` が固定環であることの実質)。
* `hres` : ★**完全分岐**(剰余体が伸びない)。落とすと偽。
* `hadj` : `O_{K′} = O_{K′′}[π′]`(原典 Lemma 5.11)。
* `[FaithfulSMul G B]` : ★落とすと偽。

★`ℕ∞` の除算・切り詰め引き算は 1 箇所も書いていない(`lean-idioms.md` #102)。 -/
theorem addVal_map_eq_card_mul {C B : Type*} [CommRing C] [IsDomain C]
    [IsDiscreteValuationRing C] [CommRing B] [IsDomain B] [IsDiscreteValuationRing B] [Algebra C B]
    {G : Type*} [Group G] [MulSemiringAction G B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H]
    {π : B} (hπ : Irreducible π)
    (hinj : Function.Injective (algebraMap C B))
    (hinvC : ∀ ρ ∈ H, ∀ x : C, ρ • algebraMap C B x = algebraMap C B x)
    (hfix : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    (hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hadj : Algebra.adjoin C ({π} : Set B) = ⊤) (c : C) :
    addVal B (algebraMap C B c) = (Nat.card H : ℕ∞) * addVal C c := by
  classical
  obtain ⟨ϖ, hϖ⟩ := exists_irreducible C
  have he0 : addVal B (algebraMap C B ϖ) ≠ 0 := fun h =>
    hϖ.not_isUnit (isUnit_of_isUnit_map hinj hinvC hfix (addVal_eq_zero_iff.mp h))
  have hne : algebraMap C B ϖ ≠ 0 := fun h => hϖ.ne_zero (hinj (by rw [h, map_zero]))
  obtain ⟨e, he⟩ : ∃ e : ℕ, addVal B (algebraMap C B ϖ) = (e : ℕ∞) := by
    obtain ⟨e, he⟩ := ENat.ne_top_iff_exists.mp (fun h => hne (addVal_eq_top_iff.mp h))
    exact ⟨e, he.symm⟩
  have hstep1 := addVal_map_eq_mul_addVal (B := B) hϖ he0
  -- `e ≤ |H|` : 共役の積(ノルム)が `C` に入り、その付値が `|H|`
  obtain ⟨ϖ₀, hϖ₀⟩ := hfix (∏ τ : H, (τ : G) • π) fun ρ hρ => smul_prod_conj π hρ
  have hval : addVal B (algebraMap C B ϖ₀) = (Nat.card H : ℕ∞) := by
    rw [hϖ₀]; exact addVal_prod_conj hπ
  have hcard : (Nat.card H : ℕ∞) = (e : ℕ∞) * addVal C ϖ₀ := by
    rw [← hval, hstep1 ϖ₀, he]
  have hcard0 : (Nat.card H : ℕ∞) ≠ 0 := by
    exact_mod_cast (Nat.card_pos (α := H)).ne'
  have hv1 : 1 ≤ addVal C ϖ₀ := by
    refine Order.one_le_iff_ne_zero.2 fun h => hcard0 ?_
    rw [hcard, h, mul_zero]
  have hle1 : e ≤ Nat.card H := by
    have : (e : ℕ∞) ≤ (Nat.card H : ℕ∞) := by
      rw [hcard]
      calc (e : ℕ∞) = (e : ℕ∞) * 1 := (mul_one _).symm
        _ ≤ (e : ℕ∞) * addVal C ϖ₀ := mul_le_mul_right hv1 _
    exact_mod_cast this
  -- `|H| ≤ e` : 完全分岐から monic 関係を作り、相異なる共役を根として数える
  have hfg : (⊤ : Submodule C B).FG := by
    have hfin := fg_adjoin_of_finite (R := C) (Set.finite_singleton π) fun x hx => by
      rw [Set.mem_singleton_iff] at hx
      subst hx
      exact isIntegral_of_fixed hfix
    rwa [hadj, Algebra.top_toSubmodule] at hfin
  obtain ⟨cc, hrel⟩ := exists_pow_eq_sum hπ ((mem_maximalIdeal ϖ).2 hϖ.not_isUnit) he hres hfg
  have hle2 : Nat.card H ≤ e :=
    card_le_of_pow_eq_sum hinvC hrel fun τ₁ τ₂ h => smul_conj_injective hadj hinvC τ₁ τ₂ h
  have hEq : (e : ℕ∞) = (Nat.card H : ℕ∞) := by
    exact_mod_cast le_antisymm hle1 hle2
  rw [hstep1 c, he, hEq]

/-- ★**ノルム `∏_{τ∈H} τ•π′` は固定環 `C` の素元である**(`v_{K′′}(ϖ₀) = 1`)。

★本体の当初の段取りは「これを先に示して `e = |H|` を出す」だったが、逆になった ——
`e = |H|` を §4(中山の補題 + 根の数え上げ)で先に出し、これはその**系**として落ちる。
`|H| = e · v_{K′′}(ϖ₀)` と `e = |H| > 0` から `v_{K′′}(ϖ₀) = 1`。 -/
theorem addVal_eq_one_of_map_eq_prod_conj {C B : Type*} [CommRing C] [IsDomain C]
    [IsDiscreteValuationRing C] [CommRing B] [IsDomain B] [IsDiscreteValuationRing B] [Algebra C B]
    {G : Type*} [Group G] [MulSemiringAction G B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H]
    {π : B} (hπ : Irreducible π)
    (hinj : Function.Injective (algebraMap C B))
    (hinvC : ∀ ρ ∈ H, ∀ x : C, ρ • algebraMap C B x = algebraMap C B x)
    (hfix : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    (hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hadj : Algebra.adjoin C ({π} : Set B) = ⊤)
    {ϖ₀ : C} (hϖ₀ : algebraMap C B ϖ₀ = ∏ τ : H, (τ : G) • π) :
    addVal C ϖ₀ = 1 := by
  have h1 : addVal B (algebraMap C B ϖ₀) = (Nat.card H : ℕ∞) := by
    rw [hϖ₀]; exact addVal_prod_conj hπ
  have h2 : (Nat.card H : ℕ∞) = (Nat.card H : ℕ∞) * addVal C ϖ₀ := by
    conv_lhs => rw [← h1]
    rw [addVal_map_eq_card_mul hπ hinj hinvC hfix hres hadj ϖ₀]
  have hne : ϖ₀ ≠ 0 := by
    intro h
    rw [h, map_zero, addVal_zero] at h1
    exact (ENat.coe_ne_top _) h1.symm
  obtain ⟨m, hm⟩ : ∃ m : ℕ, addVal C ϖ₀ = (m : ℕ∞) := by
    obtain ⟨m, hm⟩ := ENat.ne_top_iff_exists.mp fun h => hne (addVal_eq_top_iff.mp h)
    exact ⟨m, hm.symm⟩
  rw [hm, ← Nat.cast_mul] at h2
  have h3 : Nat.card H * 1 = Nat.card H * m := by
    rw [mul_one]; exact_mod_cast h2
  have h4 : (1 : ℕ) = m := Nat.eq_of_mul_eq_mul_left (Nat.card_pos (α := H)) h3
  rw [hm, ← h4, Nat.cast_one]

/-- ★★★**節点の名前どおりの形 `e(K′/K′′) = |H|`** ——
`C` の素元 `ϖ` の `B` での付値がちょうど `|H|`。 -/
theorem addVal_map_uniformizer_eq_card {C B : Type*} [CommRing C] [IsDomain C]
    [IsDiscreteValuationRing C] [CommRing B] [IsDomain B] [IsDiscreteValuationRing B] [Algebra C B]
    {G : Type*} [Group G] [MulSemiringAction G B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H]
    {π : B} (hπ : Irreducible π)
    (hinj : Function.Injective (algebraMap C B))
    (hinvC : ∀ ρ ∈ H, ∀ x : C, ρ • algebraMap C B x = algebraMap C B x)
    (hfix : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    (hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hadj : Algebra.adjoin C ({π} : Set B) = ⊤) {ϖ : C} (hϖ : Irreducible ϖ) :
    addVal B (algebraMap C B ϖ) = (Nat.card H : ℕ∞) := by
  rw [addVal_map_eq_card_mul hπ hinj hinvC hfix hres hadj ϖ, addVal_uniformizer hϖ, mul_one]

/-- `A` の像が `C` の像に含まれるなら `A[π] = ⊤ ⟹ C[π] = ⊤`。
★原典 Lemma 5.11 の `O_{K′} = O[π′]` から `O_{K′} = O_{K′′}[π′]` を出すための橋。 -/
theorem adjoin_eq_top_of_image {A C B : Type*} [CommRing A] [CommRing C] [CommRing B]
    [Algebra A B] [Algebra C B] {π : B}
    (hAC : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a)
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤) : Algebra.adjoin C ({π} : Set B) = ⊤ := by
  rw [eq_top_iff]
  rintro b -
  have hb : b ∈ Algebra.adjoin A ({π} : Set B) := by rw [hadj]; trivial
  induction hb using Algebra.adjoin_induction with
  | mem x hx => rw [Set.mem_singleton_iff.1 hx]; exact Algebra.subset_adjoin rfl
  | algebraMap a =>
      obtain ⟨c, hc⟩ := hAC a
      rw [← hc]
      exact Subalgebra.algebraMap_mem _ c
  | add x y _ _ hx hy => exact Subalgebra.add_mem _ hx hy
  | mul x y _ _ hx hy => exact Subalgebra.mul_mem _ hx hy

/-- ★★★★**Yoshida 2008 Lemma 6.8 —— `hram` を仮定に置かない形**。

Y9 `card_mul_ramIndex_eq_sum_ramIndex` の仮定 `hram` を本ファイルの主結論で供給した。
★Y9 のファイルは書き換えていない(差し替えは上流の判断)。

Y9 に対して増えた仮定は
`hπ'`(`π′` が素元)/ `hinj` / `hfixC`(`B^H ⊆ C` の像)/ `hres`(完全分岐)/
`hAC`(`O ⊆ O_{K′′}`)の 5 本で、いずれも原典 §6.2 の設定
(`K′/K` 完全分岐、`C = O_{K′′} = B^H`)が満たすものである。 -/
theorem card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing
    {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H] {π' π'' : B}
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [Algebra C B]
    [MulSemiringAction G C]
    (hcomp : ∀ (ρ : G) (c : C), algebraMap C B (ρ • c) = ρ • algebraMap C B c)
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c)
    (hπ' : Irreducible π')
    (hinj : Function.Injective (algebraMap C B))
    (hfixC : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    (hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hAC : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a)
    {ϖ : C} (hϖ : algebraMap C B ϖ = π'')
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A ({π''} : Set B))
    (σ : G) :
    (Nat.card H : ℕ∞) * ramIndex ϖ σ = ∑ τ : H, ramIndex π' (σ * (τ : G)) := by
  have hinvC : ∀ ρ ∈ H, ∀ x : C, ρ • algebraMap C B x = algebraMap C B x := fun ρ hρ x => by
    rw [← hcomp ρ x, hHtriv ρ hρ x]
  exact card_mul_ramIndex_eq_sum_ramIndex hcomp hHtriv
    (addVal_map_eq_card_mul hπ' hinj hinvC hfixC hres (adjoin_eq_top_of_image hAC hadj))
    hϖ hadj hfix σ

end ABC3.Found.PGC
