/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38KummerExists
import ABC3.Found.GenEll.Thm38Alpha
import ABC3.Found.GaloisRep.TateGaloisStab
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Meta.Claim

/-!
# 第 1160 ブロック —— **`α` が像に入る段の橋**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19–p.20。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か

`Theorem 3.8` を `Found` で閉じるのに要るものは**ちょうど 2 つ**である。

| # | 中身 | 状態 |
|---|---|---|
| I | `l`-巡回の読みの橋（`HasLCyclicVelu` ⟷ `HasLCyclicJ`） | ☆`Skeleton/GenEll/LCyclicReading.lean`（進捗枠 **2.7 / 3**、残り重み 4） |
| II | **`α` が mod `l` 像に入ること** | ☆**本ファイル**（進捗枠 **3 / 3** の部品が揃った、残りは配管 1） |

★位相の側（像が閉部分群）は第 772、群論の側（`α` と安定直線から `SL₂`）は
`§9-992` と `Lemma 3.1, (iv)` で済んでいる（`imageContainsSL2J_of_alpha'`）。

## ★★★★★★★★在庫（すでに `Found` にあるもの）

| 定理 | 内容 | 第 |
|---|---|---|
| `exists_sigma_smul_root_of_val` | `l ∤ v(q)` なら `σ(π) = η·π` なる `σ` が **`AdjoinRoot (Xˡ − C q)` の上に**ある | 994 |
| `not_lth_power_of_val` | `l ∤ v(q)` なら `q` は `l` 乗でない | 994 |
| `not_lth_power_of_not_dvd_val` | `l ∤ v_K(q)` なら `q` は `K(ζ_l)` でも `l` 乗でない（分岐指数は `l` と素） | 994 |
| `sigma_acts_as_alpha` | `σ(ζ)=ζ`・`σ(π)=ζπ` なら `σ(ζᵃπᵇ) = ζ^{a+b}πᵇ` | 993 |
| `upperM_one_mulVec` | `α = (1 1 / 0 1)` は `(a,b) ↦ (a+b, b)` | 993 |
| `tatePhi_pointMap` | Tate 一意化は `Point.map` と可換（同変性の群の形） | — |
| `tate_stab_of_pointStab` | 点の側の安定性を単数の側へ移す段 | — |

☆すなわち**両端は揃っている**——`Lˣ/qℤ` の側の `σ` と、行列 `α` の作用と、
Tate 一意化の同変性である。★足りないのは**それらを繋ぐ 3 本**である。

## ★★★★★★★★★★★★節点（進捗枠 **3 / 3** の部品が揃った）

| # | 節点 | 内容 | 重み |
|---|---|---|---|
| 1 | `tateTorsion_basis_zeta_pi` | ★**核はちょうど `lℤ × lℤ`（第 1161-1162、`zeta_pi_mem_zpowers_iff`）**。☆残るのは全射の側（`E[l]` が実際にこの 2 元で生成されること）と `Q` の無限位数の紐付け | ★★**閉じた**（第 1165 `zeta_pi_basis`） | 8 → **0** |
| 2 | `alpha_mem_localImage` | ★**座標の側（第 1164）と群論の側（第 1170 `exists_zeta_pi_of_torsion`）は揃った**。★★**第 1171 で完成**（`torsion_eq_phi_zeta_pi`） | 8 → **0** |
| 3 | `alpha_mem_globalImage` | ★**制限準同型は第 1167 で取れた**（`restrictLocalHom`・`restrictLocalHom_commutes`）。☆残るのは `L̄ ↪ M` を実際に取る段（`IsAlgClosed.lift`）と `galRep` との合成 | ★★**閉じた**（第 1168 `globalOfLocalHom`）。★★**第 1169 で完備**（`mem_map_of_mem_map_comp`） | 6 → **0** |

☆総重み 22 → **1**（第 1161-1165 で節点 1、第 1167-1169 で節点 3、
第 1170-1171 で節点 2 が閉じた）。

### ★★★★★★★★★★★★第 1171 で 3 節点がすべて揃った

`torsion_eq_phi_zeta_pi`——`Φ : Additive (G ⧸ ⟨Q⟩) ≃+ P` を Tate 一意化と読めば、
`l • Φ(c) = 0` なら `cˡ = 1`、したがって `c = [ζᵃπᵇ]`。
☆すなわち **`E[l]` の元はすべて `Φ [ζᵃπᵇ]` の形**である。

★★★**残るのは具体の `TateSetup`・`galRep` に当てはめる配管だけ**である。

### ★★★★★★★★第 1200——Kummer 拡大を Tate 設定の体として使える

`irreducible_X_pow_sub_C_of_not_pow`・`fact_irreducible_X_pow_sub_C_of_not_pow`
（`Found/GenEll/Thm38KummerExists.lean`）。

☆`l` 素数で `q` が `l` 乗でなければ `Xˡ − C q` は既約なので
`AdjoinRoot (Xˡ − C q)` は**体**になる。
★`TateSetup` は `[Field K]` を要求するので、これがないと当てはめられなかった。
☆`q` が `l` 乗でないことは `l ∤ v(q)` から出る（`not_lth_power_of_val`、第 994）。

### ★★★★★★★★★★★★★★★★第 1211——**単数の側は済んだ**

`exists_units_sigma_kummer`（`Found/GenEll/Thm38KummerExists.lean`）:

    l ∤ v(q) → ∃ π : Kˣ, ∃ σ : Kˣ →* Kˣ,
      πˡ = q ∧ σ ζ = ζ ∧ σ π = ζ · π     （K = AdjoinRoot (Xˡ − C q)）

☆これは `tate_sigma_coord_alpha`（第 1174）が**受け取る形そのもの**である
——第 994 の代数同型を `Units.map` で単数に落とし、
第 1178 の根 `π` と組み合わせた。
★`σ ζ = ζ` は `σ` が `K₀`-代数同型だから（`AlgEquiv.commutes`）。

### ★★★★★★★★★★★★★★★★★★★★第 1212——**底変換は要らなかった**

`kummer_sigma_coord_alpha`（`Found/GenEll/Thm38KummerAlpha.lean`）:

    σ (ζᵃ πᵇ) = ζ^{a'} π^{b'} Qⁿ  ⇒  l ∣ (a + b − a')  かつ  l ∣ (b − b')

☆仮説は**基礎体 `K₀` の側だけ**である——
`l ∤ v(q)`・`q` が無限位数・`ζ₀` が原始 `l` 乗根。

★★★**`TateSetup` を `K` へ底変換する段は要らなかった**
——`sigma_coord_alpha` は**任意の可換群**で成り立ち、要るのは
「`ζ` が原始 `l` 乗根」「`πˡ = Q`」「`Q` が無限位数」だけだからである。

☆すなわち `σ` は `mod l` の座標で `α = (1 1 / 0 1)` として作用する。

### ★★★★★★★★★★★★★★★★第 1227-1228——**行列への読み替えが取れた**

| 定理 | 内容 | 第 |
|---|---|---|
| `addEquiv_limTors_of_addEquiv` | ★Tate 一意化で `T_l ≅ ℤ_l²` を移す | 1227 |
| `glRed_galRep_eq_of_redVec` | ★★★座標の作用から**行列**を読む | 1228 |

☆第 1228:

    ∀ x, redVec (e (galTate σ x)) = α ·ᵥ redVec (e x)
      ⇒  glRedPadic l (galRep σ) = α

★★★**これで `α ∈ (galRep の像).map (glRedPadic l)` が言える**
——局所から大域へは第 1167-1169 が運ぶ。
### ★★★★★★★★★★★★★★★★第 1229-1232——**`ζ, π` の基底は作らなくてよい**

| 定理 | 内容 | 第 |
|---|---|---|
| `exists_conj_upperOne` | ★★自明でない幂単行列は `α` に共役 | 1229 |
| `exists_gl_lift` | ★`GL₂(𝔽_l)` の元は `GL₂(ℤ_l)` に持ち上がる | 1230 |
| `basisChange_realize` | ★任意の `GL₂(ℤ_l)` の元は基底の取り替えで実現 | 1231 |
| `exists_basis_glRed_conj` | ★★★基底を取り替えれば像は**任意の共役**に取れる | 1232 |

★★★したがって **`T_l E` の `ζ, π` 適合基底（無限の塔）を作る必要はない**
——`σ` の `mod l` の行列が**幂単かつ非自明**でありさえすれば、
基底を取り替えて `α` そのものにできる。

### ★★★★★★★★★★★★★★★★★★★★第 1233-1236——**II 側の到達点**

| 定理 | 内容 | 第 |
|---|---|---|
| `redVec_galTate` | ★行列は `mod l` 座標に `mulVec` で作用（第 1228 の逆） | 1233 |
| `glRed_unipotent_of_galTate` | ★★`σ` が `mod l` で幂単 ⇒ 行列も幂単 | 1234 |
| `glRed_ne_one_of_galTate` | ★★`σ` が非自明 ⇒ 行列も `≠ 1` | 1235 |
| `exists_basis_glRed_eq_alpha` | ★★★**ある基底で行列は `α`** | 1236 |

★★★**これで `α ∈ (galRep の像).map (glRedPadic l)` が出る**。

☆繋いだのは 4 段: 幂単性（第 1234）・非自明性（第 1235）→
`α` への共役（第 1229）→ `ℤ_l` への持ち上げ（第 1230）→
基底の取り替え（第 1231-1232）。

### ★★★★★★★★★★★★★★★★★★★★★★★★第 1237——**`halpha` の形になった**

`alpha_mem_map_of_galTate`（`Found/GenEll/AlphaMemImage.lean`）:

    σ が mod l で幂単かつ非自明
      ⇒ ∃ e₀, toGL (upper 1) ∈ (galRep E.W l e₀).range.map (glRedPadic l)

★★★**これが `imageContainsSL2J_of_alpha`（在庫）の `halpha` そのもの**である。

☆したがって II 側は **`σ` の `mod l` の幂単性と非自明性**だけに帰着した
——それを Tate 一意化（第 1212・1227）で確かめればよい。

### ★★★★★★★★★★★★第 1172-1174——仮説がすべて消えた

| 仮説 | どこから出るか | 第 |
|---|---|---|
| `hQinf`（`Q` は無限位数） | `TateSetup` の `0 < v(q)` | 1172 |
| `hμ`（`μ_l = ⟨ζ⟩`） | `IsPrimitiveRoot`（整域） | 1173 |
| `hζl`・`hζprim` | `IsPrimitiveRoot`（単数群への持ち上げ） | 1174 |

★`tate_sigma_coord_alpha` と `tate_torsion_eq_phi_zeta_pi`（第 1174）が
**Tate 設定での完成形**であり、受けているのは
「`ζ` が原始 `l` 乗根」「`πˡ = q`」「`σ(ζ) = ζ`・`σ(π) = ζπ`」だけである。
☆残るのはこれらを `L_v(ζ_l, q^{1/l})` で実際に取る段（第 994 の Kummer の `σ`）と、
`galRep` の行列に読み替える段である。

### ★★★★★★★★★★★★★★★★★★★★第 1270-1274——II 側の残りは「局所体の建設」だけ

| 定理 | 内容 | 第 |
|---|---|---|
| `galTate_unipotent_of_galPoint` | ★★`T_l E` の `mod l` 幂単性は `E[l]` で確かめれば足りる | 1270 |
| `exists_galTate_ne_of_galPoint` | ★★非自明性も同様 | 1270 |
| `point_map_galPoint` | ★埋め込みは `galPoint` と可換 | 1271 |
| `torsionMap_bijective` | ★★★代数閉体の埋め込みは `E[n]` の**全単射**（個数の勘定だけ） | 1271 |
| `galPoint_unipotent_of_map`・`exists_galPoint_ne_of_map` | ★2 条件が埋め込みで降りる | 1271 |
| `tate_unipotent_of_sigma`・`tate_exists_ne_of_sigma` | ★★★**`σ` は `E[l]` に幂単・非自明**（抽象な `Φ`・`τ` で） | 1272 |
| `tate_point_unipotent`・`tate_point_exists_ne` | ★★点の側で読んだ形（同変性は `tatePhi_pointMap`） | 1273 |
| `exists_algEquiv_sigma_kummer` | ★Kummer の `σ` を**体自己同型**として取る | 1274 |

★★☆★★☆**訂正（第 1275・1277）**——下の表の B1（局所体の建設）は**要らない**。
また在庫の `tatePhi_pointMap` は `σA : K →ₐ[R] K` が恒等射しかないため
（第 1275 で証明）、点の側の同変性は `tatePhi_map`（`σR`・`σK` の対）から
書き直す必要がある。☆第 1276（`tateSigmaAct`）がその受け皿である。

### ★★★★★★★★★★★★★★★★★★★★訂正後の道筋（第 1277）

☆**`π` は体を拡げなくても作れる**——`exists_pow_eq_of_torsion_not_mu`（第 1277）:

> `l`-捩れの類が `μ_l` の像に収まらなければ、その類の代表 `x` は `x^l = Q^k`（`l ∤ k`）を満たし、
> Bezout で `π ≔ x^b Q^a` が `π^l = Q` を与える。

★★★したがって道筋は「**大域の数体を `L′ ≔ L(ζ_l, E[l], √d)` に取り替え、
その完備化 `P.adicCompletion L′` で `mkTateSetup`（在庫）を使う**」となる。

| # | 訂正後の残り | 在庫 |
|---|---|---|
| 1 | `L′` を取る（`ζ_l`・`E[l]`・捻り `d` を添加。`L/L′` は有限 Galois） | `exists_finite_point_descent`（第 1207）等 |
| 2 | `P.adicCompletion L′` で `TateSetup` と `Φ` | `mkTateSetup`・`dvrTatePhiAddEquiv`（**無条件**） |
| 3 | `E ⊗ K′ ≅ tateCurveAt q`（変数変換は `R′` の上） | `exists_tate_model`（在庫） |
| 4 | `π` を捩れから作る | **第 1277** |
| 5 | `σ` を分解群から取り、`σπ = ζπ` にする | `exists_algEquiv_sigma_kummer`（第 1274）の類推 |
| 6 | `τ = 実際の Galois 作用` | `tatePhi_map`＋`tateCurveAt_map`（在庫） |
| 7 | 大域へ運び `T_l E` に上げる | 第 1271・第 1270 |

★★**局所体の建設（完備 DVR の Eisenstein 拡大）は道から消えた**。

★★★**数学の中身は尽きた**。残っているのは次の 3 本の「局所体の建設」である:

| # | 残る節点 | 中身 | 見立て |
|---|---|---|---|
| B1 | `TateSetup` を `K ≔ L_v(ζ_l, q^{1/l})` に建てる | ★完備 DVR の**完全分岐（Eisenstein）拡大**が再び完備 DVR であること。☆第 1013-1018 は**不分岐 2 次**の場合（`f̄` 既約）で、`X^l − q` は `f̄ = X^l` なので当たらない | 既知数学。mathlib の Henselian ＋ Dedekind から組める |
| B2 | `hσv`（付値が `σ` で不変） | ☆完備体の付値の拡大は一意。★`tateDvrVal` を整閉包から作れば `τ(R′) = R′` で済む | B1 に付随 |
| B3 | Tate 曲線と悪い素点の `E` を結ぶ | ☆`E ⊗ L_v` は変数変換で `tateCurveAt q` に写る（`tateParamR`・`integralModel` は在庫） | 在庫の組み替え |

☆B1-B3 が済めば、`restrictLocalHom`（第 1167）と第 1271 で大域へ運び、
第 1270 で `T_l E` に上げ、`alpha_mem_map_of_galTate`（第 1237）に入る。

### ★★★★★★★★★★★★★★★★★★★★第 1276-1282——局所の `σ` は完全に手に入った

| 定理 | 内容 | 第 |
|---|---|---|
| `tateSigmaAct` | ★`τ ≔ Φ ∘ σU の誘導写像 ∘ Φ⁻¹`（同変性は定義から） | 1276 |
| `tateSigmaAct_unipotent`・`_exists_ne` | ★★★`τ` は `E[l]` で幂単かつ非自明 | 1276 |
| `exists_pow_eq_of_torsion_not_mu` | ★`μ_l` からはみ出す類があれば `Q` は `l` 乗 | 1277 |
| `exists_torsion_not_mu` | ★`l`-捩れが `l` 個より多ければはみ出す | 1278 |
| `exists_pi_of_phi` | ★★★**`E[l]` が載っていれば `π` は作れる**（体を拡げない） | 1279 |
| `not_mem_range_of_not_dvd_vAdd` | ★`l ∤ v(Q)` なら `π` は基礎体に無い | 1280 |
| `exists_algEquiv_move` | ★基礎体に無い元はどれかの `σ` が動かす | 1280 |
| `sigma_pi_eq_zeta_mul` | ★★`σ(π) = ζπ` で `ζ` は原始 `l` 乗根 | 1281 |
| `exists_sigma_zeta_pi` | ★★★**`σ`・`ζ`・`π` の三つ組**（第 1276 の入力が完成） | 1282 |

☆すなわち **`l ∤ v(Q)` と「基礎体が `μ_l` を含む」と「`E[l]` が載っている」の 3 つから、
`τ` の幂単性・非自明性が出る**。★どれも局所体の建設を要しない。

### ★★★★★★★★残り 2 本（2026-09-02 時点）

| # | 残る節点 | 中身 | 道具 |
|---|---|---|---|
| 6 | `tateSigmaAct` ＝ 実際の Galois 作用 | ☆`tatePhi_map`（`σR`・`σK` の対）を点の言葉に直す。曲線は `σ` で固定される部分環 `R₀` の上に置く | `tateCurveAt_map`・`exists_point_image_eq`（`subst` の技） |
| 8 | 悪い素点での組み立て | ☆`L′ ≔ L(ζ_l, E[l], √d)` を取り、`P.adicCompletion L′` で `mkTateSetup`・`dvrTatePhiAddEquiv`・`exists_tate_model` | `not_dvd_vAdd_tateParam_of_not_dvd_jExp`（在庫）が `l ∤ v(Q)` を与える |

★（7 大域への輸送は第 1271 で済み、`T_l E` への持ち上げは第 1270 で済み。）

### ★★★★★★★★★★★★★★★★★★★★★★★★第 1283-1290——鎖は「局所データを与える」まで詰んだ

| 定理 | 内容 | 第 |
|---|---|---|
| `pointCoords_tatePhi_sigma` | ★★★同変性の**正しい形**（座標の言葉、`σR`・`σK` の対） | 1283 |
| `map_vcX_fixed`・`map_vcY_fixed` | ★座標は写像と可換 | 1284 |
| `galAct`・`pointCoords_galAct` | ★★`σ` が固定する曲線の上の点の作用 | 1285 |
| `tateSigmaAct_eq_galAct` | ★★★**`tateSigmaAct` ＝ 実際の Galois 作用** | 1286 |
| `exists_galPoint_conditions_global` | ★★★局所の 2 条件は大域へ降りる | 1287 |
| `galAct_vcPoint` | ★★`σ` の作用は変数変換と可換 | 1288 |
| `galAct_unipotent_ne_of_tate` | ★★★**局所側の完成**（Tate モデルで幂単・非自明） | 1289 |
| `unipotent_ne_of_variableChange` | ★★2 条件は `E ⊗ K` に戻る | 1290 |

☆**鎖はこう繋がる**:

    第 1289（Tate モデルで幂単・非自明）
      → 第 1290（`E ⊗ K` へ）
      → 第 1287（大域 `L̄` へ）
      → 第 1270（`T_l E` へ）
      → 第 1237（`α` が mod `l` 像に入る）
      → 第 1249（`ImageContainsSL2J`）

### ★★★★★★★★残るのは節点 8（局所データを与える段）だけ

☆第 1289 が受け取るのは次だけであり、どれも「作れる」ことが分かっている:

| # | 与えるもの | 出どころ |
|---|---|---|
| a | `L′ ≔ L(ζ_l, E[l], √d)`（有限 Galois）と素点 `P`、`K ≔ P.adicCompletion L′` | `exists_finite_subextension`（第 1195）・`numberField_and_towers`（第 1222） |
| b | `E ⊗ K` の分裂乗法還元 → `q`・変数変換 `C` | `exists_tate_model`・`tateParamR`（在庫） |
| c | `TateSetup`・`Φ` | `mkTateSetup`・`dvrTatePhiAddEquiv`（在庫、**無条件**） |
| d | `π`（`π^l = Q`） | 第 1279（`E[l]` が `K` に載っていればよい） |
| e | `l ∤ v(Q)` | `not_dvd_vAdd_tateParam_of_not_dvd_jExp`（在庫）＋ `PrimeToLocalHeights` |
| f | `σ`（`σζ = ζ`・`σπ = ζπ`）と `σR`・`hσv` | 第 1282（`IsGalois K₀ K` が要る） |

★★`ζ_l` は**大域の体に添加してよい**——`ImageContainsSL2J` は部分群の像で言えば十分だからである。

### ★★☆★★☆第 1291——残るのは (f) の 1 命題だけ

☆(a)-(e) はすべて在庫の型に落ちる。★**残るのは (f)**:

> `K ≔ P.adicCompletion L′` は `K₀ ≔ p.adicCompletion L` の上で**正規**である
> （同値: `K` は `K₀` と `L′` の合成体である）

☆これは数論の標準的な事実だが **mathlib にも本プロジェクトにも無い**（2026-09-02 実測）。
★これが無いと「`π` を動かす `K` の自己同型」が取れない
——`K₀(π)/K₀` だけなら分裂体なので正規だが、`TateSetup` を `K₀(π)` に建てるには
完備 DVR の Eisenstein 拡大が要り、第 1277 で消したはずの道に戻ってしまう。

☆(e) について念のため: `l ∤ v(Q)` は**基礎体 `K₀` の付値**についてでよい
（第 1280・1282 が使うのは `v₀`）ので、`L′/L` の分岐指数が `l` で割れても構わない。
### ★★★★★★★★★★★★★★★★★★★★★★★★第 1292-1297——道が短くなった（節点 (f) は消えた）

☆第 1291 で「完備化の正規性が要る」と測ったが、**幂単性はもっと安く出る**ことが分かった:

| 定理 | 内容 | 第 |
|---|---|---|
| `sq_sub_one_eq_zero_of_det_one_of_fixed` | ★★★`det = 1` ＋ 固定ベクトル ⇒ `(M−1)² = 0` | 1292 |
| `galTate_unipotent_of_matrix` | ★★行列が幂単 ⇒ `h2` | 1293 |
| `exists_fixed_vec_of_galPoint_eq` | ★★固定点 ⇒ 固定ベクトル | 1294 |
| `det_glRed_eq_one_of_fixed_root` | ★★`ζ_l` を固定 ⇒ `det ≡ 1` | 1295 |
| `galTate_unipotent_of_fixed_root_point` | ★★★**`h2` は「`ζ_l` と固定点」だけから出る** | 1296 |
| `addOrderOf_tatePhi_zeta` | ★★`Φ[ζ]` は位数ちょうど `l`（`μ_l ⊂ E(K₀)`） | 1297 |

☆**新しい道**（Tate 一意化の同変性も局所体の拡大も要らない）:

| 条件 | 出どころ |
|---|---|
| `h2`（幂単） | `ζ_l ∈ L`（大域で添加）＋ `μ_l ⊂ E(K₀)`（第 1297、**基礎局所体の上の Tate 一意化だけ**） |
| `h1`（非自明） | `E[l] ⊄ E(K₀)`（第 1279＋`l ∤ v(Q)`）＋ `exists_algEquiv_move`（第 1280） |

★★★**第 1291 で名指しした節点 (f)（`P.adicCompletion L′` の正規性）は道から消えた**。
☆残るのは配管だけである:

* (P1) `Φ[ζ]` を Tate モデルから `E ⊗ L̄` へ運ぶ（第 1288・1290・1271 の道具）
* (P2) `E[l] ⊄ E(K₀)` から `σ` を取る（第 1280・1279・`torsion_card`）
* (P3) 悪い素点で `TateSetup` を建てる（`mkTateSetup`・`exists_tate_model`、在庫）
* (P4) `ζ_l ∈ L` にするための大域の底変換（`ImageContainsSL2J` は部分群の像で十分）

### ★★★★★★★★★★★★★★★★★★★★★★★★第 1298-1303——II 側は 3 つの入力だけに帰着した

| 定理 | 内容 | 第 |
|---|---|---|
| `exists_galPoint_ne_of_coord_not_mem` | ★座標が基礎体に無ければ動かされる | 1298 |
| `galPoint_rhPoint_eq` | ★基礎体の点は `σ` に固定される | 1299 |
| `galTate_unipotent_of_rational` | ★★★基礎体に `ζ_l` と位数 `l` の点があれば**どの `σ` も**幂単 | 1300 |
| `exists_galTate_ne_of_coord_not_mem` | ★★`h1` は座標が基礎体に無い点から | 1301 |
| `exists_galPoint_fixed_of_map` | ★★固定点も埋め込みで降りる | 1302 |
| `galTate_h2_h1_of_fixed_moved` | ★★★**`h2`・`h1` の一括** | 1303 |

☆**II 側が要求するのはこの 3 つだけ**になった:

| # | 入力 | 出どころ |
|---|---|---|
| 1 | `σ` が原始 `l` 乗根 `ζ` を固定 | `ζ_l ∈ L`（大域で添加してよい——`ImageContainsSL2J` は部分群の像で十分） |
| 2 | `σ` が固定する位数 `l` の点 | `μ_l ⊂ E(K₀)`（第 1297）＋ 第 1302 で大域へ |
| 3 | `σ` が動かす `l`-捩れ点 | `E[l] ⊄ E(K₀)`（`l ∤ v(Q)`、第 1279・1298）＋ 第 1271 で大域へ |

★★★**どれも「基礎局所体 `K₀ = p.adicCompletion L` の上の Tate 一意化」から出る**
——`mkTateSetup`・`dvrTatePhiAddEquiv`・`exists_tate_model` はすべて在庫（無条件）である。
☆第 1291 で名指しした「完備化の正規性」も、第 1274 の「Eisenstein 拡大」も、
**どちらも道から消えた**。

### ★★★★★★★★★★★★★★★★★★★★★★★★第 1304-1311——鎖は「Tate 一意化を建てる」まで詰んだ

| 定理 | 内容 | 第 |
|---|---|---|
| `card_torsion_le_of_not_dvd` | ★★`l ∤ v(Q)` ⇒ 基礎体の `l`-捩れは `l` 個以下 | 1304 |
| `exists_galPoint_ne_of_not_rational` | ★基礎体の点でなければ動かされる | 1305 |
| `exists_torsion_not_rational` | ★★個数の勘定（`l² > l`） | 1306 |
| `exists_sigma_h2_h1_of_local` | ★★★局所データから `h2`・`h1` | 1307 |
| `exists_h2_h1_global_of_local` | ★★★局所の `σ` から大域へ（輸送） | 1308 |
| `exists_local_fixed_moved` | ★★局所で固定点・動く点を同時に | 1309 |
| `baseChange_baseChange`・`galPoint_pointEquivOfEq` | ★塔の突き合わせ | 1310 |
| `exists_h2_h1_global_of_localData` | ★★★★**橋の完成形** | 1311 |

☆**鎖はこう繋がる**:

    悪い素点で `TateSetup`（`mkTateSetup`・`dvrTatePhiAddEquiv`、在庫・無条件）
      → `μ_l ⊂ E(L_v)`（第 1297）と `l`-捩れ ≤ `l` 個（第 1304）
      → 第 1311（局所データ → 大域の `h2`・`h1`）
      → 第 1237（`α` が mod `l` 像に入る）
      → 第 1249（`ImageContainsSL2J`）

### ★★★★★★★★残るのは 2 つ（2026-09-02 時点）

| # | 残り | 中身 |
|---|---|---|
| a | 悪い素点で `TateSetup` を建て、`μ_l` と `l`-捩れの個数を `E` の言葉にする | `exists_tate_model`（在庫）で Tate モデルへ移し、`vcPoint` で戻す |
| b | `ζ_l ∈ L` にするための大域の底変換 | `ImageContainsSL2J` は部分群の像で十分（制限準同型） |

★★**未知の数学は 1 本も無い**——どちらも在庫の組み替えである。

### ★★★★★★★★★★★★★★★★★★★★★★★★★★★★第 1312-1314——鎖は 1 本に繋がった

    `mkTateSetup`・`dvrTatePhiAddEquiv`（在庫、**無条件**）
      → 第 1314（局所の 2 入力: 位数 `l` の点・`l`-捩れ ≤ `l` 個）
      → 第 1313（変数変換で `E ⊗ L_v` へ）
      → 第 1311／1312（大域の `h2`・`h1`）
      → 第 1237（`α` が mod `l` 像に入る）
      → 第 1249（`ImageContainsSL2J`）

☆**第 1312 で `ζ_l ∈ L`（大域の底変換）が要らなくなった**——
必要なのは `σ` が `ζ` を固定することだけで、それは `ζ ∈ L_v` から出る。
★`L_v` を「`L(ζ_l)` の素点での完備化」に取ればよく、
`restrictLocalHom` は `L` の上で働くので `SSCurve` を底変換しなくてよい。

### ★★★★★★★★残るのは「悪い素点での組み立て」だけ

☆具体的には、`E` が `p` で乗法還元（`HasMultRed`）かつ `l ∤ v_p(q)`
（`PrimeToLocalHeights`）のとき:

1. `L(ζ_l)` の `p` 上の素点 `P` を取り、`L_v ≔ P.adicCompletion` とする（`ζ_l ∈ L_v`）
2. `E ⊗ L_v` の分裂乗法還元（必要なら不分岐 2 次捻り——第 1007-1029 の在庫）
3. `exists_tate_model`（在庫）で `q` と変数変換 `C` を取る
4. `mkTateSetup`・`dvrTatePhiAddEquiv`（在庫、無条件）で `S`・`Φ`
5. `not_dvd_vAdd_tateParam_of_not_dvd_jExp`（在庫）で `l ∤ v(Q)`

★★**未知の数学は 1 本も無い**——すべて在庫の組み替えである。

### ★★★★★★★★★★★★★★★★★★★★★★★★★★★★第 1315-1317——局所側は「分裂乗法還元」まで詰んだ

| 定理 | 内容 | 第 |
|---|---|---|
| `exists_vc_tateModel` | ★`E` と Tate モデルを結ぶ変数変換（`K` の上で） | 1315 |
| `exists_point_order_of_vc′`・`card_le_of_vc′` | ★曲線の等式つきで移す | 1316 |
| `local_inputs_of_split` | ★★★★**分裂乗法還元から局所の 2 入力** | 1317 |

☆**鎖の現状**:

    分裂乗法還元（悪い素点）＋ 極小性 ＋ `ζ ∈ K` ＋ `l ∤ v(Q)`
      → 第 1317（`W ⊗ L_v` についての 2 入力）
      → 第 1311／1312（大域の `h2`・`h1`）
      → 第 1237（`α` が mod `l` 像に入る）
      → 第 1249（`ImageContainsSL2J`）

### ★★★★★★★★残るのは 3 つの「在庫の呼び出し」

| # | 残り | 在庫 |
|---|---|---|
| i | 悪い素点で分裂乗法還元を用意する（必要なら不分岐 2 次捻り） | 第 1007-1029 の `UnramQuad` 一式 |
| ii | `l ∤ v(Q)` を `PrimeToLocalHeights` から出す | `not_dvd_vAdd_tateParam_of_not_dvd_jExp` |
| iii | `ζ ∈ K` にするため素点を `L(ζ_l)` の上で取る | `adicCompletion`（在庫）＋ 第 1312（底変換は不要） |

★★**未知の数学は 1 本も無い**。

### ★★★★★★★★第 1170 で節点 2 の群論の側が閉じた

`exists_zeta_pi_of_torsion`——`G ⧸ ⟨Q⟩` の `l`-捻れはすべて `[ζᵃπᵇ]` である。
☆Tate 一意化 `Φ : G ⧸ ⟨Q⟩ ≃+ E.Point` は加法同型なので、
`E[l]` の元は `Φ [ζᵃπᵇ]` の形に書ける。
★★**残るのは `Φ` を介して `E[l]` の言葉に直す配管だけ**である。
★★**残るのは具体の `TateSetup`・`galRep` への当てはめだけ**である。

### ★★★★★★★★第 1169 で節点 3 が完備した

`mem_map_of_mem_map_comp`——`f` を `galRep`、`g` を `globalOfLocalHom`、
`r` を `glRedPadic l` と読むと、これが
`imageContainsSL2J_of_alpha'` の `halpha` を**局所から受け取る形**である。

### ★★★★★★★★第 1168 で節点 3 が閉じた

`globalOfLocalHom`——`M` が代数閉体なら `L̄ = AlgebraicClosure L` は
`IsAlgClosed.lift` で `M` に埋め込めるので、その埋め込みで
`restrictLocalHom` を当てればよい。
★★`Gal(M/L_v) →* Gal(L̄/L)` ができた——**局所で実現される行列は大域の像にも入る**。

### ★★★★★★★★第 1167 で取れたもの（節点 3 の核）

`Found/GenEll/Thm38Decomposition.lean` の **`restrictLocalHom`**——
塔 `L → L_v → M` と `L` 上正規な中間体 `E` に対し

    `σ : M ≃ₐ[L_v] M`  →  `σ|_L : M ≃ₐ[L] M`  →  `σ|_E : E ≃ₐ[L] E`

の準同型。☆`restrictLocalHom_commutes` で埋め込みと可換であることも取った。

★★**単射性は要らない**——必要なのは「局所の元の像が大域の像に**含まれる**」
ことだけで、それは準同型があれば出る。

### ★★★★★★★★第 1161 で取れたもの

`Found/GenEll/Thm38ZetaPi.lean` の **`zeta_pi_indep`**——
`ζᵃ·πᵇ ∈ ⟨Q⟩` なら `l ∣ a` かつ `l ∣ b` である（★無条件）。

☆機構: `ζᵃπᵇ = Qⁿ` を `l` 乗すると `Qᵇ = Q^{ln}`、`Q` は無限位数なので `b = ln`。
すると `πᵇ = (πˡ)ⁿ = Qⁿ` なので `ζᵃ = 1`、`ζ` の位数は `l` なので `l ∣ a`。
★`orderOf_eq_of_primitive`（原始 `l` 乗根の位数はちょうど `l`）も一緒に取った。

★第 1162 で逆向き（`l ∣ a` かつ `l ∣ b` なら `ζᵃπᵇ ∈ ⟨Q⟩`）も取り、
**核がちょうど `lℤ × lℤ`** であることが iff の形で出た（`zeta_pi_mem_zpowers_iff`）。
☆これが「`([ζ], [π])` は `E[l]` の `ℤ/l`-基底である」の正確な形であり、
節点 2 が消費する界面である。

★第 1163 でさらに `zeta_pi_coord_eq_iff`——
`[ζ^{a₁}π^{b₁}] = [ζ^{a₂}π^{b₂}] ⟺ a₁ ≡ a₂, b₁ ≡ b₂ (mod l)`——を取った。
☆これで座標は `(ZMod l)²` の元として**一意に**決まる。
★`sigma_acts_as_alpha`（第 993）の `(a,b) ↦ (a+b, b)` と合わせれば、
`σ` の行列がちょうど `α = (1 1 / 0 1)` になる。

★★第 1164 でそれを組んだ `sigma_coord_alpha` を取った——
`σ(ζᵃπᵇ)` を座標 `(a', b')` で書いたなら必ず `a' ≡ a + b`、`b' ≡ b`（mod `l`）。
☆**これで座標の側は全部揃った**。

★★★第 1165 で**全射の側**（`zeta_pi_span`: `xˡ ∈ ⟨Q⟩` なら `x = ζᵃπᵇ`）も取り、
`zeta_pi_basis`（生成＋一意性）で**節点 1 が閉じた**。
☆残るのは Tate 一意化で `E[l]` を `Lˣ/qℤ` の `l`-捻れと同定する段（節点 2 の残り）と、
分解群経由で大域へ移す段（節点 3）だけである。

★節点 1-3 が閉じれば `imageContainsSL2J_of_alpha'` の `halpha` が埋まる。

### ★★★★節点 3 の道（測定）

`Gal(L̄_v/L_v) → Gal(L̄/L)` は**単射ではない**が、`L̄ ↪ L̄_v` を取れば
制限写像 `Gal(L̄_v/L_v) → Gal(L̄/L)` が定義され、その像が分解群である。
★`galRep` は `Gal(L̄/L)` の表現なので、局所の `σ` の像は大域の像に含まれる。
☆mathlib の `AlgHom.restrictNormal` / `IsAlgClosed.lift` が素材である。

### ☆節点 1 が本体である理由

`E[l]` の `ℤ/l`-基底を `(ζ, π)` に取ることが、`α = (1 1 / 0 1)` という**具体的な行列**を
出す唯一の道である。★`Lˣ/qℤ` の `l`-捩れは `μ_l · π^ℤ / qℤ` であり、
`πˡ = q` なので `π` の類は位数 `l`、`ζ` の類も位数 `l`、両者は独立である。

## ★★★★★★何が `Theorem 3.8` に残るか（第 1160 の総括）

    I（`l`-巡回の読み、残り 4） ＋ II（本ファイル、残り 1） ＝ **5**

☆`Theorem 3.8` 以外（`Corollary 4.3`・`4.4`）は本定理の上にしか立たない。
★`Theorem 2.1`（§2）は曲線の Riemann–Roch が mathlib に無いので別筋である。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta

/-! ## ★出典の紐付け(`.src`) -/

def alphaBridgeFrame.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(α が mod l 像に入る段の橋——3 節点の枠)",
    sectionId := "genell-thm-3-8" }

def alphaBridgeFrame.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_sigma_smul_root_of_val(σ の存在、第 994、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_sigma_smul_root_of_val") 1,
    .citation "[ABC3]" "sigma_acts_as_alpha(行列表示、第 993、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.sigma_acts_as_alpha") 1,
    .citation "[ABC3]" "imageContainsSL2J_of_alpha'(残る仮説は α だけ、第 772、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.imageContainsSL2J_of_alpha'") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1160）の測定**——`Theorem 3.8` を `Found` で閉じるのに" ++
       "要るのはちょうど 2 つ: (I) `l`-巡回の読みの橋（`LCyclicReading`、残り重み 5）と " ++
       "(II) `α` が mod `l` 像に入ること（本ファイル、第 1161 で 22 → 17）。" ++
       "☆両端——`Lˣ/qℤ` の側の `σ`（第 994）と行列 `α` の作用（第 993）と " ++
       "Tate 一意化の同変性（`tatePhi_pointMap`）——はすべて揃っている。" ++
       "★足りないのは (1) `E[l]` の基底を `(ζ, π)` に取る段、" ++
       "(2) 局所の mod `l` 像が `α` を含むこと、" ++
       "(3) 分解群経由で大域の像へ移す段の 3 本である。") 22 ]

end ABC3.Skeleton.GenEll
